import os
import ROOT
import pandas as pd
import seaborn as sns
import matplotlib

path = "./analysis/UL2018_v2/mutau_eff/DY/"

version = ["2p1","2p5"]
fakes   = ["jet","mu", "e"]
wp      = ["L","M","T"]
source  = ["R","E", "J","M"]

        
dic_2p1 = {"tot": [0,0,0,0]}
dic_2p5 = {"tot": [0,0,0,0]}

for w in wp:
    for f in fakes:
        dic_2p1["idVS"+f+"_"+w] = [0,0,0,0]
        dic_2p5["idVS"+f+"_"+w] = [0,0,0,0] 

df_2p1 =pd.DataFrame(dic_2p1)
df_2p5 =pd.DataFrame(dic_2p5)

df_2p1.index = source
df_2p5.index = source

dict_name = {}
print(os.listdir(path)[0])
f = os.path.join(path, os.listdir(path)[0])
if os.path.isfile(f):
   _file0 = ROOT.TFile.Open(f)
   h = _file0.Get("cutflow")
   for i in range(3,79):
       name = h.GetXaxis().GetBinLabel(i)
       dict_name[name] = i

print(dict_name)
#phr = "tau_idVSmu_L_2p1_R"
#words = phr.split('_')

#key  = words[1]+"_"+words[2]
#vers = words[3]
#src  = words[4]

#print(vers+src)
#print(key)


for filename in os.listdir(path):

    f = os.path.join(path, filename)
    # checking if it is a filif os.path.isfile(f):
    if os.path.isfile(f):
	print(f)
	_file0 = ROOT.TFile.Open(f)
	t = _file0.Get("cutflow")

        for s in source: 
            df_2p1.loc[s,"tot"] += t.GetBinContent(dict_name["tau_"+s])
            df_2p5.loc[s,"tot"] += t.GetBinContent(dict_name["tau_"+s])
            for w in wp:
                for f in fakes:
                    content = t.GetBinContent(dict_name["tau_idVS"+f+"_"+w+"_2p1_"+s])
                    print("-------------------------------")
                    print("tau_idVS"+f+"_"+w+"_2p1_"+s)
                    print("v2p1: " +str( content))
                    df_2p1.loc[s,"idVS"+f+"_"+w] += content

                    content = t.GetBinContent(dict_name["tau_idVS"+f+"_"+w+"_2p5_"+s])
                    print("v2p5: " + str(content))
                    df_2p5.loc[s,"idVS"+f+"_"+w] += content
print(df_2p1)
print(df_2p5)

#do ratios
for col in df_2p1.columns:
    if col != "tot": 
       df_2p1.loc[:,col] = df_2p1[col]/df_2p1["tot"]
       df_2p5.loc[:,col] = df_2p5[col]/df_2p5["tot"]

print(df_2p1)

#plot table and save
df_2p1.pop("tot")
df_2p5.pop("tot")
norm = matplotlib.colors.Normalize(-1,1)
colors = [[norm(-1.0), "white"],
          [norm( 1.0), "white"]]
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
fig = sns.heatmap(df_2p1, annot=True, cbar=False).get_figure()
fig.savefig("efficiencies_v2p1.png")

fig2 = sns.heatmap(df_2p5, annot=True, cbar=False).get_figure()
fig2.savefig("efficiencies_v2p5.png")

