"""
Prepare data to show the results in tables
"""

import sys
import numpy as np
import pandas as pd
import os

BASE = sys.argv[1] if len(sys.argv)>=2 else 'rfc'
M = int(sys.argv[2]) if len(sys.argv)>=3 else 100

EXP_PATH = "../out/"
DATASETS = {
            'haberman':'Haberman',
            'breast-cancer':'Breast-C',
            'tictactoe':'TicTacToe',
            'bank-notes':'Banknote',
            'kr-vs-kp':'kr-vs-kp',
            'spam':'Spambase',
            'mushroom':'Mushroom',
            'adults':'Adult',
            }

PREC = 3

PRETTY_MAP = {
    "best_mv_risk":"$L(h_{best})$",
    "unf_mv_risk":"$L(\\MV_{u})$",
    "lam_mv_risk":"$L(\\MV_{\\rho_\\lambda})$",
    "tnd_mv_risk":"$L(\\MV_{\\rho_{\\TND}})$",
    "cctnd_mv_risk":"$L(\\MV_{\\rho_{\\CCTND}})$",
    "ccpbb_mv_risk":"$L(\\MV_{\\rho_{\\CCPBB}})$",
    "lam_PBkl":"$\\FO(\\rho_\\lambda)$",
    "tnd_TND":"$\\TND(\\rho_{\\TND})$",
    "cctnd_CCTND":"$\\CCTND(\\rho_{\\CCTND})$",
    "ccpbb_CCPBB":"$\\CCPBB(\\rho_{\\CCPBB})$",
    "lam":"$\\FO$",
    "tnd":"$\\TND$",
    "cctnd":"$\\CCTND$",
    "ccpbb":"$\\CCPBB$",
    "gibbs":"$\\E_\\rho[L]$",
    "tandem":"$\\E_{\\rho^2}[L]$",
    "bmu":"$\mu$",
}

def comparison_table(tp='risk'):
    path = "table/"+base+"/optimize/"
    out_fname = path+tp+"_table.tex"
    if not os.path.isdir(path):
        os.makedirs(path)

    opts = {
        ("risk","rfc"):       ["unf","lam","tnd","cctnd","ccpbb"],
        ("risk","mce"):       ["unf","best","lam","tnd","cctnd","ccpbb"],
        ("bound","rfc"):      [("lam","PBkl"),("tnd","TND"),("cctnd","CCTND"),("ccpbb","CCPBB")],
        ("bound","mce"):      [("lam","PBkl"),("tnd","TND"),("cctnd","CCTND"),("ccpbb","CCPBB")],
    }[(tp,base)]
    if tp=='risk':
        opts = [(o,"mv_risk") for o in opts]
    
    copts = [pre+"_"+suf for pre,suf in opts]
    
    if hl1=="all":
        hl1 = copts

    with open(out_fname, 'w') as fout:
        # Header
        fout.write("\\begin{tabular}{l"+"c"*len(opts)+"}\\toprule\n")
        fout.write("Data set")
        for i,col in enumerate(copts):
            fout.write(" & "+PRETTY_MAP[col])
        fout.write(" \\\\\n")
        fout.write("\\midrule\n")

        for ds in DATASETS:
            if (base == 'rfc' and ds == 'Protein'):
                continue
            df = pd.read_csv(EXP_PATH+ds+"-"+str(M)+"-bootstrap-iRProp.csv",sep=";")
            df_mean = df.mean()
            df_std  = df.std()
            
            # Highlight indices
            v1 = np.min(df_mean[hl1]) if len(hl1)>0 else -1
            v2 = np.min(df_mean[hl2]) if len(hl2)>0 else -1
            v1 = str(round(v1,PREC))
            v2 = str(round(v2,PREC))

            fout.write("\\dataset{"+RENAME.get(ds,ds)+"}")
            for i,col in enumerate(copts):
                fval = df_mean[col]
                val = str(round(fval,PREC))
                std = str(round(df_std[col],PREC))
                s = val + " ("+std+")"
                if col in hl1 and val==v1:
                    s = "\\textbf{"+s+"}"
                if col in hl2 and val==v2:
                    s = "\\underline{"+s+"}"
                fout.write(" & "+s)
            fout.write(" \\\\\n")

        fout.write("\\bottomrule\n") 
        fout.write("\\end{tabular}\n")
    
#optimized_comparison_table('risk', base=BASE, hl2=["lam_mv_risk","tnd_mv_risk","cctnd_mv_risk","ccpbb_mv_risk"])
#optimized_comparison_table('bound', base=BASE, hl2=["tnd_TND","cctnd_CCTND","ccpbb_CCPBB"])


"""
Making tables comparing FW and FWEL bounds:
"""
def FW_FWEL_table():
    path = "table/"
    out_fname = path+"FW_FWEL_table.tex"
    if not os.path.isdir(path):
        os.makedirs(path)
    
    opts = ["PBkl","MGG","Skl"]
    rep_opts = {"bound":"Bound","Lntest":"Risk","bestSigma2":"Var"}
    
    with open(out_fname, 'w') as fout:
        # Header
        fout.write("\\begin{tabular}{|l|l|"+"ll|"*len(opts)+"}\n")
        fout.write("\\hline\n")
        fout.write("\\multirow{2}{*}{Dataset} & ")
        for o in opts:
            fout.write(" & \\multicolumn{2}{c|}{$\\"+o+"$}")
        fout.write(" \\\\\n")
        fout.write(" &"+" & \\multicolumn{1}{c}{$\\FW$} & \\multicolumn{1}{c|}{$\\FWEL$}"*len(opts)+" \\\\\n")
        fout.write("\\hline\n")
        

        for ds in DATASETS.keys():
            dfFW = pd.read_csv(EXP_PATH+ds+"-"+"FW.csv",sep=",")
            dfFW_mean = dfFW.mean()
            dfFWEL = pd.read_csv(EXP_PATH+ds+"-"+"FWEL.csv",sep=",")
            dfFWEL_mean = dfFWEL.mean()
            
            fout.write("\\multirow{3}{*}{\\dataset{"+DATASETS.get(ds,ds)+"}}")
            for r in rep_opts.keys():
                fout.write("&"+rep_opts[r])
                for o in opts:
                    cell = r+"_"+o
                    fout.write("&"+str(round(dfFW_mean[cell],PREC))+"&"+str(round(dfFWEL_mean[cell],PREC)))
                fout.write(" \\\\\n")
            fout.write("\\hline\n")
        fout.write("\\end{tabular}\n")

FW_FWEL_table()


"""
Making tables about values of FWEL:
"""
def values(method="FWEL"):
    path = "table/"
    out_fname = path+method+"_values.tex"
    if not os.path.isdir(path):
        os.makedirs(path)
    
    Left = ["ExL1","RefL1"]
    Right = ["Term1","RefTerm1"]
    opts = ["PBkl","MGG","Skl"]
    
    with open(out_fname, 'w') as fout:
        # Header
        fout.write("\\begin{tabular}{|l|l|"+"c|"*(len(Left)+1)+"}\n")
        fout.write("\\hline\n")
        fout.write("Dataset & Methods")
        fout.write(" & $\E_{\\rho} [\\Delta_{\\hat L}(h,h_{S_1},S_2)]$ / $\\B(\E_{\\rho} [\\Delta_L(h,h_{S_1})])$ ")
        fout.write(" & $\\hat L(h_{S_1},S_2)$ / $\B(L(h_{S_1}))$ ")
        fout.write(" & $\E_{\\rho} [\\Delta^+_{\\hat L}(h,h_{S_1},S_2)]$ / $\E_{\\rho} [\\Delta^-_{\\hat L}(h,h_{S_1},S_2)]$")
        fout.write(" \\\\\n")
        fout.write("\\hline\n")
        
        for ds in DATASETS.keys():
            dfFWEL = pd.read_csv(EXP_PATH+ds+"-"+"FWEL.csv",sep=",")
            dfFWEL_mean = dfFWEL.mean()
            
            fout.write("\\multirow{3}{*}{\\dataset{"+DATASETS.get(ds,ds)+"}}")
            for o in opts:
                fout.write("&"+o)
                for i in range(len(Left)):
                    fout.write("&"+str(round(dfFWEL_mean[Left[i]+"_"+o],PREC))+" / "+str(round(dfFWEL_mean[Right[i]+"_"+o],PREC)))
                fout.write("&")
                if o=="Skl":
                    fout.write(str(round(dfFWEL_mean["ExL1P"],PREC))+" / "+str(round(dfFWEL_mean["ExL1M"],PREC)))
                fout.write(" \\\\\n")
            fout.write("\\hline\n")
        fout.write("\\end{tabular}\n")

values("FWEL")

# NOT USED
def TND_Ben_comparison_table(base='rfc'):
    path = "table/"+base+"/optimize/"
    if not os.path.isdir(path):
        os.makedirs(path)
    
    prec = 5
    opts = ["tnd", "mu", "bern"]
    cols = ["dataset", "c", "d"]
    for opt in opts:
        if opt == "tnd":
            cols += [opt+suf for suf in ["_KL", "_gibbs", "_tandem", "_tnd", "_TandemUB"]]
        elif opt == "mu":
            cols += [opt+suf for suf in ["_KL", "_gibbs", "_tandem", "_MU", "_muTandemUB", "_bmu"]]
        elif opt == "bern":
            cols += [opt+suf for suf in ["_KL", "_gibbs", "_tandem", "_bern", '_mutandem_risk', '_vartandem_risk', "_varUB", "_bernTandemUB", "_bmu", "_bg", "_bl"]]
    rows = []
    for ds in DATASETS:
        if (base == 'rfc' and ds == 'Protein'):
            continue
        df = pd.read_csv(EXP_PATH+ds+"-"+str(M)+"-bootstrap-iRProp.csv",sep=";")
        df_mean = df.mean()
        df_std  = df.std()
        
        row = [ds, df_mean["c"], df_mean["d"]]
        for opt in opts:
            if opt == "tnd":
                row += [df_mean[opt+suf] for suf in ["_KL", "_gibbs", "_tandem", "_tnd", "_TandemUB"]]
            elif opt == "mu":
                row += [df_mean[opt+suf] for suf in ["_KL", "_gibbs", "_tandem", "_MU", "_muTandemUB", "_bmu"]]
            elif opt == "bern":
                row += [df_mean[opt+suf] for suf in ["_KL", "_gibbs", "_tandem", "_bern", '_mutandem_risk', '_vartandem_risk', "_varUB", "_bernTandemUB", "_bmu", "_bg", "_bl"]]            
        rows.append(row)
    
    pd.DataFrame(data=rows, columns=cols).round(prec).to_csv(path+"mu_comparison.csv", sep=",", index=False)

#TND_Ben_comparison_table(base=BASE)
