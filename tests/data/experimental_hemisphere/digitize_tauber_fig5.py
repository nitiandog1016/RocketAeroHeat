"""人工选取原图独立圆点；仅生成草稿数据，不采样关联式曲线。"""
import csv,json,hashlib,pathlib,math
from PIL import Image,ImageDraw
p=pathlib.Path(__file__).parent
points=[(154,69),(283,67),(377,211),(386,233),(446,243),(454,310),
        (497,285),(542,359),(547,322),(559,304),(566,340),(571,412),
        (645,595),(665,514),(684,504),(709,464),(710,632),(768,699),
        (784,561),(819,664),(823,751),(847,801),(813,822),(814,838)]
# 原页200%渲染后裁剪(630,125,1085,650)，再放大2倍。
# 以下坐标属于裁剪放大图；横轴线性0--90度，纵轴log10(.1)--log10(1)。
x0,x90,y1,y01=66.,881.,66.,958.
rows=[]
for i,(x,y) in enumerate(points,1):
    phi=(x-x0)*90./(x90-x0)
    ratio=10.**(-(y-y1)/(y01-y1))
    rows.append(dict(point_id=i,crop_pixel_x=x,crop_pixel_y=y,
        original_page_pixel_x=630+x/2.,original_page_pixel_y=125+y/2.,
        body_angle_deg=phi,q_over_q_stagnation=ratio,
        angle_extraction_halfwidth_deg=4.*90./(x90-x0),
        q_ratio_extraction_relative_halfwidth=10.**(4./(y01-y1))-1.,
        symbol='circle',series_identity='not_verified_from_original_legend',
        status='digitization_draft_conditions_incomplete'))
with (p/'tauber_fig5_circle_points_draft.csv').open('w',newline='',encoding='utf-8') as f:
    w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
im=Image.open(p/'tauber_fig5_crop.png').convert('RGB'); d=ImageDraw.Draw(im)
for i,(x,y) in enumerate(points,1):
    d.ellipse((x-6,y-6,x+6,y+6),outline=(220,35,35),width=2)
    d.text((x+7,y-6),str(i),fill=(170,25,25))
im.save(p/'tauber_fig5_selected_points.png')
metadata=dict(status='draft',source='NASA TM-101077 (1989), printed p8, Fig5, PDF page12',
    official_pdf_url='https://ntrs.nasa.gov/api/citations/19890017745/downloads/19890017745.pdf',
    official_pdf_sha256=hashlib.sha256((p/'19890017745.pdf').read_bytes()).hexdigest(),
    underlying_original='Kemp, Rose, Detra (1959), JAS26(7),421-430, DOI10.2514/8.8128, original Fig3',
    original_source_access='Publisher403; secondary mirror OCR read, original experiment legend and wallT not independently verified',
    point_count=len(rows),point_type='visually selected experimental circles, not theoretical curve samples',
    axes=dict(x='s/R in degrees, 0 to90 linear',y='q/q_stagnation, .1 to1 log10',
              crop_calibration=dict(x0=x0,x90=x90,y1=y1,y01=y01)),
    extraction_uncertainty='2 source page pixels (4 enlarged crop pixels), only nominal click/calibration uncertainty',
    experimental_uncertainty='Not reconstituted. Original reported gauge covers5--10deg; experiments aggregate multiple Ms/p1.',
    eligibility_for_independent_experiment_comparison=False,
    missing=['verified original legend and run identity per point','wall temperature','shock-tube gas upstreamT',
             'local inviscid conditions or reconstructed shock-tube upstream flow','catalysis and chemistry compatibility',
             'gauge averaging and experimental uncertainty'])
(p/'tauber_fig5_points_metadata.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2),encoding='utf-8')
assert len(rows)==24 and all(0<r['body_angle_deg']<90 and 0.1<r['q_over_q_stagnation']<=1 for r in rows)
print('PASS:24 selected circles, linear/log axes decoded, retained raw pixel provenance; draft only.')
