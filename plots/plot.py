"""
Prepare data to show the results in tables
"""

import sys
import numpy as np
import pandas as pd
import os

EXP_PATH = "../out/"
DATASETS = {
            'haberman':'Haberman',
            'breast-cancer':'Breast-C',
            'tictactoe':'TicTacToe',
            'bank-notes':'Banknote',
            'kr-vs-kp':'kr-vs-kp',
            'spam':'Spambase',
            'svmguide1':'Svmguide1',
            'mushroom':'Mushroom',
            'adults':'Adult',
            }

PREC = 3

""" Plot error and bounds for several data sets """
def multi_bounds():
    path = "figure/datasets/"
    if not os.path.isdir(path):
        os.makedirs(path)

    for ds in DATASETS.keys():
        df = pd.read_csv(EXP_PATH+ds+"-"+"Avg-Ex-IP.csv",sep=",")
        df_mean, df_std = df.mean(), df.std()
        bounds = ["PBkl", "PBUB", "PBSkl"]
        with open(path+ds+".tex", "w") as f:
            for i, bnd in enumerate(bounds):
                f.write("\\addplot["+bnd+", Bound]coordinates {("+str(i+1)+","+str(df_mean["bound_"+bnd])+") +- (0,"+str(df_std["bound_"+bnd])+")};\n")
            for i, bnd in enumerate(bounds):
                f.write("\\addplot["+bnd+", MyRisk]coordinates {("+str(i+1)+","+str(df_mean["Lntest_"+bnd])+") +- (0,"+str(df_std["Lntest_"+bnd])+")};\n")

multi_bounds()

"""
Making tables comparing Avg bounds:

def All_table(rep_opt="bound"):
    # rep_opt: "bound", "Lntest", "bestSigma2"
    path = "table/"
    out_fname = path+"All_"+rep_opt+"_table.tex"
    if not os.path.isdir(path):
        os.makedirs(path)
    
    with open(out_fname, 'w') as fout:
        # Header
        if rep_opt == "Lntest":
            NOC = 8
        else:
            NOC = 7
        fout.write("\\begin{tabular}{|l|"+"c|"*NOC+"}\n")

        fout.write("\\hline\n")
        fout.write(" Dataset ")
        if rep_opt == "Lntest":
            fout.write("& $h_S$ ")
        fout.write("& $\\PBkl$ & $\\MGG$ & $\\PBkl_{\\FW}$ & $\\MGG_{\\FW}$ & $\\PBkl_{\\Avg}$ & $\\MGG_{\\Avg}$ & $\\Skl_{\\Avg}$ \\\\\n")
        fout.write("\\hline\n")
        
        for ds in DATASETS.keys():
            df = pd.read_csv(EXP_PATH+ds+"-"+".csv",sep=",")
            dfFW = pd.read_csv(EXP_PATH+ds+"-"+"FW.csv",sep=",")
            dfAvg = pd.read_csv(EXP_PATH+ds+"-"+"Avg.csv",sep=",")
            df_mean = df.mean()
            dfFW_mean = dfFW.mean()
            dfAvg_mean = dfAvg.mean()
            
            rep_means = [df_mean[rep_opt+"_PBkl"], df_mean[rep_opt+"_MGG"], dfFW_mean[rep_opt+"_PBkl"], dfFW_mean[rep_opt+"_MGG"]\
                            , dfAvg_mean[rep_opt+"_PBkl"], dfAvg_mean[rep_opt+"_MGG"], dfAvg_mean[rep_opt+"_Skl"]]
            # Highlight indices
            v1 = np.min(rep_means)
            v2 = np.min(rep_means[-3:])
            v1 = str(round(v1,PREC))
            v2 = str(round(v2,PREC))

            fout.write("\\dataset{"+DATASETS.get(ds,ds)+"}")
            if rep_opt == "Lntest":
                fout.write("& "+str(round(df_mean["LnERMtest"],PREC)))
            for i in range(len(rep_means)):
                val = str(round(rep_means[i],PREC))
                s = val
                if val==v1:
                    s = "\\textbf{"+s+"}"
                if i > (NOC-3) and val==v2:
                    s = "\\underline{"+s+"}"
                fout.write(" & "+s)
            fout.write(" \\\\\n")
        fout.write("\\hline\n")
        fout.write("\\end{tabular}\n")
"""
#All_table("bound")
#All_table("Lntest")
#All_table("bestSigma2")


"""
Making tables comparing Vanilla and FW of PBkl and MGG bounds:

def Vanilla_FW_table():
    path = "table/"
    out_fname = path+"Vanilla_FW_table.tex"
    if not os.path.isdir(path):
        os.makedirs(path)
    
    cols = {"bound":"bound", "bestSigma2":"$\\sigma^*$", "Lntrain":"$\\E_{\\rho}[\\hat L(h,S)]$", "KL":"$\\KL$"}
    opts = ["PBkl","MGG"]
    opts2 = ["Vanilla", "FW"]
    
    with open(out_fname, 'w') as fout:
        # Header
        fout.write("\\begin{tabular}{|l|l|"+"cccc|"*len(opts2)+"}\n")
        fout.write("\\hline\n")
        fout.write("\\multirow{2}{*}{Dataset} & ")
        for o in opts2:
            fout.write(" & \\multicolumn{"+str(len(cols.keys()))+"}{c|}{"+o+"}")
        fout.write(" \\\\\n")
        fout.write(" &")
        for _ in range(len(opts2)):
            for c in cols.keys():
                fout.write(" & \\multicolumn{1}{c|}{"+cols[c]+"}")
        fout.write(" \\\\\n")
        fout.write("\\hline\n")
        
        for ds in DATASETS.keys():
            df = pd.read_csv(EXP_PATH+ds+"-"+".csv",sep=",")
            df_mean = df.mean()
            dfFW = pd.read_csv(EXP_PATH+ds+"-"+"FW.csv",sep=",")
            dfFW_mean = dfFW.mean()
            
            fout.write("\\multirow{"+str(len(opts))+"}{*}{\\dataset{"+DATASETS.get(ds,ds)+"}}")
            for o in opts:
                fout.write("&"+o)
                for c in cols.keys():
                    fout.write("&"+str(round(df_mean[c+"_"+o],PREC)))
                for c in cols.keys():
                    fout.write("&"+str(round(dfFW_mean[c+"_"+o],PREC)))
                fout.write(" \\\\\n")
            fout.write("\\hline\n")
        fout.write("\\end{tabular}\n")
"""
#Vanilla_FW_table()

"""
Making tables comparing FW and FWEL bounds:

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
"""
#FW_FWEL_table()


"""
Making tables about values of FWEL:

def FWEL_values():
    path = "table/"
    out_fname = path+"FWEL_values.tex"
    if not os.path.isdir(path):
        os.makedirs(path)
    
    Left = ["ExL1","RefL1"]
    Right = ["ExTerm1","RefTerm1"]
    opts = ["PBkl","MGG","Skl"]
    
    with open(out_fname, 'w') as fout:
        # Header
        fout.write("\\begin{tabular}{|l|l|l|"+"c|"*(len(Left)+1)+"}\n")
        fout.write("\\hline\n")
        fout.write("Dataset & Methods")
        fout.write(" & $\\Delta_{\\ell}=0$")
        fout.write(" & $\E_{\\rho} [\\Delta_{\\hat L}(h,h_{S_1},S_2)]$ / $\\B(\E_{\\rho} [\\Delta_L(h,h_{S_1})])$ ")
        fout.write(" & $\\hat L(h_{S_1},S_2)$ / $\\B(L(h_{S_1}))$ ")
        fout.write(" & $\E_{\\rho} [\\Delta^+_{\\hat L}(h,h_{S_1},S_2)]$ / $\E_{\\rho} [\\Delta^-_{\\hat L}(h,h_{S_1},S_2)]$")
        fout.write(" \\\\\n")
        fout.write("\\hline\n")
        
        for ds in DATASETS.keys():
            df = pd.read_csv(EXP_PATH+ds+"-FWEL.csv",sep=",")
            df_mean = df.mean()
            
            fout.write("\\multirow{3}{*}{\\dataset{"+DATASETS.get(ds,ds)+"}}")
            for o in opts:
                fout.write("& $\\"+o+"$")
                fout.write("&"+str(round(df_mean["ExL1_0rate_"+o]*100,1))+" \\%")
                for i in range(len(Left)):
                    fout.write("&"+str(round(df_mean[Left[i]+"_"+o],PREC))+" / "+str(round(df_mean[Right[i]+"_"+o],PREC)))
                fout.write("&")
                if o=="Skl":
                    fout.write(str(round(df_mean["ExL1P"],PREC))+" / "+str(round(df_mean["ExL1M"],PREC)))
                fout.write(" \\\\\n")
            fout.write("\\hline\n")
        fout.write("\\end{tabular}\n")
"""
#FWEL_values()

