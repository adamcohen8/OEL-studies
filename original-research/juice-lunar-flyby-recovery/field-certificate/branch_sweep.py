from branch_cover import *
a=Arc("juice")
index=int(sys.argv[1])
configs=json.loads((HERE/"global.json").read_text())
row=configs[index]
results=[]
for N,threshold in [(12,7000.),(12,10000.),(12,14000.),(12,20000.),(32,9000.),(32,12000.)]:
    r=solve(a,row,N,threshold,physical=index>=2)
    results.append(r);save("sweep-"+str(index)+".json",results)
    print(index,N,threshold,r["diameter_deg"],r["hit_radius_deg"],r["far_radius_deg"],flush=True)
