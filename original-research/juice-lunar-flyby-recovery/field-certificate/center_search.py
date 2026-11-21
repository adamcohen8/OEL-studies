from branch_cover import *
a=Arc("juice");index=int(sys.argv[1]);row=json.loads((HERE/"global.json").read_text())[index]
out=[]
for shift in (.5,1.,1.5,2.):
    for radius in (12000.,17000.):
        q=solve(a,row,24,radius,directional=True,shift=shift)
        out.append(dict(shift=shift,**q));save("centers-"+str(index)+".json",out)
        print(index,shift,radius,q["diameter_deg"],flush=True)
