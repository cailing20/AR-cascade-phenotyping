"""One-page vector/PDF study wheel with radial labels and source-qualified points."""
from pathlib import Path
import argparse, collections, csv, json, math, html
from openpyxl import load_workbook
import pypdfium2
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor
from reportlab.pdfbase.pdfmetrics import stringWidth

REPO = Path(__file__).resolve().parent
DEFAULT_WORKBOOK = REPO / "Supplementary_Table_Figure_4_carrier_evidence.xlsx"
O = None

W,H,CX,CY=1500,1640,750,735
SCALE=540/W
def xy(a,r):a=math.radians(a);return CX+r*math.cos(a),CY+r*math.sin(a)
COLORS={'Positive':'#C46245','Protective':'#26877C','Null':'#FFFFFF','Mixed':'#D4AC40'}
class Art:
    def __init__(self):
        self.svg=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" role="img" aria-label="Carrier phenotype evidence by domain, gene and measurement category"><rect width="{W}" height="{H}" fill="white"/>']
        self.c=canvas.Canvas(str(O/'carrier_study_wheel.pdf'),pagesize=(W*SCALE,H*SCALE),invariant=1)
        self.c.setTitle('Carrier phenotyping across disease domains');self.c.scale(SCALE,SCALE)
    def line(self,x1,y1,x2,y2,color='#DCE4E8',width=1):
        self.svg.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" y2="{y2:.3f}" stroke="{color}" stroke-width="{width}"/>')
        self.c.setStrokeColor(HexColor(color));self.c.setLineWidth(width);self.c.setDash();self.c.line(x1,H-y1,x2,H-y2)
    def circle(self,x,y,r,fill=None,stroke='#DCE4E8',width=1,dash=False):
        self.svg.append(f'<circle cx="{x:.3f}" cy="{y:.3f}" r="{r:.3f}" fill="{fill or "none"}" stroke="{stroke or "none"}" stroke-width="{width}"'+(' stroke-dasharray="3 2"' if dash else '')+'/>')
        self.c.setFillColor(HexColor(fill or '#FFFFFF'));self.c.setStrokeColor(HexColor(stroke or fill or '#FFFFFF'));self.c.setLineWidth(width);self.c.setDash(3,2) if dash else self.c.setDash();self.c.circle(x,H-y,r,stroke=bool(stroke),fill=bool(fill))
    def poly(self,pts,fill,stroke=None,width=1,dash=False):
        self.svg.append('<polygon points="'+' '.join(f'{x:.3f},{y:.3f}' for x,y in pts)+f'" fill="{fill}" stroke="{stroke or "none"}" stroke-width="{width}"'+(' stroke-dasharray="3 2"' if dash else '')+'/>')
        p=self.c.beginPath();p.moveTo(pts[0][0],H-pts[0][1])
        for x,y in pts[1:]:p.lineTo(x,H-y)
        p.close();self.c.setFillColor(HexColor(fill));self.c.setStrokeColor(HexColor(stroke or fill));self.c.setLineWidth(width);self.c.setDash(3,2) if dash else self.c.setDash();self.c.drawPath(p,stroke=bool(stroke),fill=1)
    def sector(self,a,z,inner,outer,fill):
        ts=[a+(z-a)*i/80 for i in range(81)];self.poly([xy(t,inner) for t in ts]+[xy(t,outer) for t in reversed(ts)],fill)
    def text(self,x,y,s,size=22,color='#243F4E',anchor='start',rotate=0,bold=False,italic=False):
        font='Helvetica-BoldOblique' if bold and italic else 'Helvetica-Bold' if bold else 'Helvetica-Oblique' if italic else 'Helvetica'
        self.svg.append(f'<text x="{x:.3f}" y="{y:.3f}" font-family="Arial, sans-serif" font-size="{size}" fill="{color}" text-anchor="{anchor}"'+(' font-weight="bold"' if bold else '')+(' font-style="italic"' if italic else '')+(f' transform="rotate({rotate:.3f} {x:.3f} {y:.3f})"' if rotate else '')+'>'+html.escape(s)+'</text>')
        c=self.c;c.saveState();c.translate(x,H-y);c.rotate(-rotate);c.setFillColor(HexColor(color));c.setFont(font,size)
        (c.drawCentredString if anchor=='middle' else c.drawRightString if anchor=='end' else c.drawString)(0,0,s);c.restoreState()
    def finish(self):
        self.svg.append('</svg>');(O/'carrier_study_wheel.svg').write_text(''.join(self.svg));self.c.showPage();self.c.save()

def _bool(value, row, name):
    if str(value).strip().lower() in {"yes", "true", "1"}:
        return True
    if str(value).strip().lower() in {"no", "false", "0"}:
        return False
    raise ValueError(f"Figure marks row {row}: {name} must be Yes or No")


def _count(value, row, name):
    if value is None or value == "":
        return None
    try:
        number = int(value)
        if float(value) != number or number < 1:
            raise ValueError
        return number
    except (TypeError, ValueError):
        raise ValueError(f"Figure marks row {row}: {name} must be a positive integer or blank") from None


def read_marks(workbook, config):
    """Read the explicit plotted marks, not cohort totals from Evidence details."""
    book = load_workbook(workbook, read_only=True, data_only=True)
    if "Figure marks" not in book.sheetnames:
        raise ValueError("Workbook lacks the 'Figure marks' plotting sheet.")
    sheet = book["Figure marks"]
    rows = sheet.iter_rows(values_only=True)
    headers = next(rows)
    required = {"point_id", "gene", "domain_fan", "canonical_study_id", "reference_id", "cohort", "tested_variants", "measurement_category", "evidence_ids", "outcome_rows", "direction", "carrier_n_min", "carrier_n_max", "abstract_only", "qualified_context", "phenotypes"}
    missing = required - set(headers)
    if missing:
        raise ValueError(f"Figure marks is missing columns: {', '.join(sorted(missing))}")
    domain_by_gene = {gene: domain for domain, genes, _, _ in config["groups"] for gene in genes}
    points, seen = [], set()
    for n, values in enumerate(rows, start=2):
        r = dict(zip(headers, values))
        if not any(v is not None and v != "" for v in values):
            continue
        for key in ("point_id", "gene", "domain_fan", "canonical_study_id", "reference_id", "measurement_category", "direction"):
            if r.get(key) is None or str(r[key]).strip() == "":
                raise ValueError(f"Figure marks row {n}: {key} is required")
            r[key] = str(r[key]).strip()
        if r["point_id"] in seen:
            raise ValueError(f"Duplicate point_id: {r['point_id']}")
        seen.add(r["point_id"])
        if r["gene"] not in domain_by_gene:
            raise ValueError(f"Figure marks row {n}: add gene {r['gene']} to config/figure.json first")
        if domain_by_gene[r["gene"]] != r["domain_fan"]:
            raise ValueError(f"Figure marks row {n}: domain for {r['gene']} disagrees with config/figure.json")
        if r["measurement_category"] not in config["categories"]:
            raise ValueError(f"Figure marks row {n}: unrecognized measurement category")
        if r["direction"] not in COLORS:
            raise ValueError(f"Figure marks row {n}: direction must be Positive, Null, Protective, or Mixed")
        for k in ("cohort", "tested_variants", "evidence_ids", "phenotypes"):
            r[k] = str(r.get(k) or "")
        r["carrier_n_max"] = _count(r.get("carrier_n_max"), n, "carrier_n_max")
        r["carrier_n_min"] = _count(r.get("carrier_n_min"), n, "carrier_n_min")
        if r["carrier_n_min"] and r["carrier_n_max"] and r["carrier_n_min"] > r["carrier_n_max"]:
            raise ValueError(f"Figure marks row {n}: carrier_n_min exceeds carrier_n_max")
        r["abstract_only"] = _bool(r["abstract_only"], n, "abstract_only")
        r["qualified_context"] = _bool(r["qualified_context"], n, "qualified_context")
        count = r["carrier_n_max"]
        bins = config["size_bins"]
        if count is None:
            r["carrier_n_bin"], r["marker_radius"] = "Unknown", 6
        else:
            idx = 0 if count < 10 else 1 if count < 100 else 2 if count < 1000 else 3 if count < 10000 else 4
            r["carrier_n_bin"], r["marker_radius"] = bins[idx]
        points.append(r)
    if not points:
        raise ValueError("The Figure marks sheet contains no points.")
    return points


def main():
    global O
    parser = argparse.ArgumentParser(description="Generate the carrier evidence wheel from the supplementary workbook.")
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK, help="Workbook with a Figure marks sheet (default: bundled detailed evidence workbook)")
    parser.add_argument("--output-dir", type=Path, default=REPO / "outputs", help="Directory for PDF, SVG, PNG, and mark summary")
    args = parser.parse_args()
    if not args.workbook.is_file():
        parser.error(f"Workbook not found: {args.workbook}")
    O = args.output_dir.resolve()
    O.mkdir(parents=True, exist_ok=True)
    config = json.loads((REPO / "config" / "figure.json").read_text())
    cats, radii, bins = config["categories"], config["radii"], config["size_bins"]
    points = read_marks(args.workbook, config)
    plotted_genes = {p["gene"] for p in points}
    groups = [[domain, [g for g in gs if g in plotted_genes], color, fill] for domain, gs, color, fill in config["groups"]]
    groups = [group for group in groups if group[1]]
    genes = [g for _, gs, _, _ in groups for g in gs]
    qualified={p['gene'] for p in points if p['qualified_context']}
    available=360-24-3*len(groups);needed={g:6.5 for g in genes}
    for gene in genes:
        for ci,cat in enumerate(cats):
            cell=[p for p in points if p['gene']==gene and p['measurement_category']==cat]
            if not cell:continue
            big=max(p['marker_radius'] for p in cell);sep=2*big+3;rows=max(1,math.floor((64-2*big)/sep)+1);cols=math.ceil(len(cell)/rows)
            needed[gene]=max(needed[gene],(cols-1)*sep/(.8*radii[ci]*math.pi/180)+3)
    assert sum(needed.values())<available,('Insufficient angular space',needed)
    pitches={g:needed[g]+(available-sum(needed.values()))/len(genes) for g in genes}
    cursor=-78;angles={};sectors=[]
    a=Art()
    for domain,gs,color,fill in groups:
        start=cursor;stop=start+sum(pitches[g] for g in gs);a.sector(start,stop,110,460,fill)
        gc=start
        for g in gs:angles[g]=gc+pitches[g]/2;gc+=pitches[g]
        mid=(start+stop)/2;rot=(mid+90)%360
        if 90<rot<270:rot-=180
        labels=domain.split(' ') if domain in ['Drug response','Psychiatry context'] else [domain]
        for j,label in enumerate(labels):
            x,y=xy(mid,657+(j-(len(labels)-1)/2)*29);a.text(x,y,label,23,color,'middle',rot,True)
        sectors.append(dict(domain=domain,genes=gs,start_angle=start,stop_angle=stop,color=color));cursor=stop+3
    for radius in [110,180,250,320,390,460]:a.circle(CX,CY,radius,None,'#D7E1E6',1)
    for gene,angle in angles.items():
        a.line(*xy(angle,110),*xy(angle,478),'#CDD9DF',1)
        x,y=xy(angle,493);norm=angle%360;left=90<norm<270
        a.text(x,y,gene+('\u2020' if gene in qualified else ''),24,anchor='end' if left else 'start',rotate=angle+180 if left else angle,bold=True,italic=True)
    # Pack complete cells, avoiding the central-first trap with five large marks.
    placed=[]
    for gene in genes:
      for ci,cat in enumerate(cats):
        cell=sorted([p for p in points if p['gene']==gene and p['measurement_category']==cat],key=lambda p:(-p['marker_radius'],p['point_id']))
        if not cell:continue
        radius=radii[ci];big=max(p['marker_radius'] for p in cell);sep=2*big+3.4;maxoff=35-big-3;n=len(cell)
        if n==1:slots=[(0,0)]
        elif n==2:slots=[(0,-sep/2),(0,sep/2)]
        else:
            cols=math.ceil(n/2);xs=[(j-(cols-1)/2)*sep*radius/(radius-maxoff) for j in range(cols)];order=sorted(range(cols),key=lambda j:(abs(xs[j]),j))
            slots=[(xs[j],dr) for j in order for dr in [-maxoff,maxoff]][:n]
        found=False
        for flipx,flipy,shift in [(fx,fy,s) for s in [0,-.08*pitches[gene],.08*pitches[gene]] for fx in [1,-1] for fy in [1,-1]]:
            trial=[]
            for p,(dx,dr) in zip(cell,slots):
                da=math.degrees(math.atan2(dx*flipx,radius))+shift;dr*=flipy;x,y=xy(angles[gene]+da,radius+dr)
                trial.append(dict(p,x=x,y=y,radial_offset=dr,angular_offset=da))
            if all(math.hypot(p['x']-q['x'],p['y']-q['y'])>=p['marker_radius']+q['marker_radius']+2.1 for p in trial for q in placed) and all(math.hypot(p['x']-q['x'],p['y']-q['y'])>=p['marker_radius']+q['marker_radius']+2.1 for j,p in enumerate(trial) for q in trial[j+1:]):
                for p,q in zip(cell,trial):p.update(q)
                placed.extend(cell);found=True;break
        if not found:raise AssertionError(('Cell packing failed',gene,cat,n))
    for p in points:
        x,y,r=p['x'],p['y'],p['marker_radius'];fill=COLORS[p['direction']]
        title=f"{p['gene']}; {p['measurement_category']}; {p['reference_id']}; carrier N {p['carrier_n_max'] or 'unknown'}; {p['direction']}; {p['phenotypes']}"
        a.svg.append(f'<g class="study-point" tabindex="0" role="button" data-point="{p["point_id"]}" aria-label="{html.escape(title,quote=True)}"><title>{html.escape(title)}</title>')
        if p['carrier_n_bin']=='Unknown':a.poly([(x,y-r),(x+r,y),(x,y+r),(x-r,y)],fill,'#334F60',1.7,p['abstract_only'])
        else:a.circle(x,y,r,fill,'#334F60',1.7,p['abstract_only'])
        a.svg.append('</g>')
    for i,r in enumerate(radii,1):
        x,y=xy(-90,r);a.circle(x,y,15,'#FFFFFF','#AEBFC8',1.2);a.text(x,y+7,str(i),22,anchor='middle',bold=True)
    pubs=len({p['canonical_study_id'] for p in points})
    a.circle(CX,CY,106,'#FFFFFF','#BDCCD5',1.6)
    a.text(52,1454,'Measurement rings',23,bold=True)
    for i,c in enumerate(cats):a.text(52,1490+i*32,f'{i+1}  {c}',22)
    a.text(525,1454,'Carrier N (categorical)',23,bold=True)
    for i,(label,r) in enumerate(bins):
        x=542+(i%2)*230;y=1502+(i//2)*51;a.circle(x,y,r,'#CFDBE1','#334F60',1.5);a.text(x+25,y+8,label,22)
    a.poly([(772,1595),(779,1602),(772,1609),(765,1602)],'#FFFFFF','#334F60',1.5);a.text(797,1610,'Unknown N',22)
    a.text(1040,1454,'Finding and source',23,bold=True)
    for i,(label,fill) in enumerate([('Difference reported',COLORS['Positive']),('Null comparison','#FFFFFF'),('Favorable / protective',COLORS['Protective']),('Mixed endpoints',COLORS['Mixed'])]):
        y=1490+i*32;a.circle(1050,y-7,7,fill,'#334F60',1.5);a.text(1070,y,label,22)
    a.circle(1050,1611,7,'#FFFFFF','#334F60',1.5,True);a.text(1070,1619,'Dashed: source limited',22)
    a.finish()
    # Use PDFium for a portable 3000-pixel PNG without a system Poppler dependency.
    pdf = pypdfium2.PdfDocument(str(O / "carrier_study_wheel.pdf"))
    page = pdf[0]
    scale = 3000 / max(page.get_width(), page.get_height())
    page.render(scale=scale).to_pil().save(O / "carrier_study_wheel.png")
    summary = {"genes": len(genes), "publications": pubs, "marks": len(points), "workbook": str(args.workbook.resolve())}
    (O / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))

if __name__=='__main__':main()
