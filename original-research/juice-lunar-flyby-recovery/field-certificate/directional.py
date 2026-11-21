from branch_cover import *
a=Arc("juice")
row=json.loads((HERE/"global.json").read_text())[int(sys.argv[1])]
out=[]
for threshold in (8000.,12000.,20000.,30000.):
    r=solve(a,row,24,threshold,directional=True)
    out.append(r);save("directional-"+sys.argv[1]+".json",out)
    print(threshold,r["diameter_deg"],r["turn_groups"],r["total_groups"],flush=True)
