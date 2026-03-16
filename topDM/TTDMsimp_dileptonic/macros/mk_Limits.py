import os
from brazil_band import getLimits, plotUpperLimits

os.environ["PYTHONNOUSERSITE"] = "1"

# Mass points
#mPhi_values = ['10','50','100','150', '200', '250', '300', '350', '400', '500', '600' '700', '800', '1000'] 
mPhi_values = ['800'] 

# Directories
datacards_dir = "datacards_files"
root_dir = "root_files"
limit_dir = "limit_files"

def run_commands(mphi):

    os.makedirs(datacards_dir, exist_ok=True)
    os.makedirs(root_dir, exist_ok=True)
    os.makedirs(limit_dir, exist_ok=True)

    comb_card_filename = os.path.join(datacards_dir, f"ttDM_CombCard_mphi_{mphi}.txt")
    root_filename = os.path.join(root_dir, f"ttDM_mphi_{mphi}.root")
    limit_filename = os.path.join(limit_dir, f"Limit_mphi_{mphi}.txt")

    combine_cmd = (
        f"combineCards.py " \
#        f"SR2b=ttdm_sr_2l_2b/phi1/datacard_TTDMsimpSpin0_s_mphi-{mphi}.txt "
        f"SR1b=ttdm_sr_2l_1b/phi1/datacard_TTDMsimpSpin0_s_mphi-{mphi}.txt " \
        f"DYCR=dycr_1b/events/datacard_TTDMsimpSpin0_s_mphi-{mphi}.txt " \
        f"ttZCR=ttZcr_inclusive/events/datacard_TTDMsimpSpin0_s_mphi-{mphi}.txt " \
        f"> {comb_card_filename}"
    )

    text2workspace_cmd = (
        f"text2workspace.py "
        f"-P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel "
        f"--PO verbose "
        f"--PO 'map=.*/TTDM*:r[1,0,10]' "
        f"{comb_card_filename} "
        f"-o {root_filename}"
    )

    combine_limit_cmd = (
        f"combine -M AsymptoticLimits "
        f"-t -1 --expectSignal 0 --cminDefaultMinimizerStrategy 0 "
        f"{root_filename} &> {limit_filename}"
    )

    os.system(combine_cmd)
    os.system(text2workspace_cmd)
    os.system(combine_limit_cmd)


def main():

    labels = []
    values = []

    for mphi in mPhi_values:

        print(f"Running mphi = {mphi}")

        run_commands(mphi)

        limit_file = os.path.join(limit_dir, f"Limit_mphi_{mphi}.txt")

        labels.append(limit_file)
        values.append(int(mphi))

    if labels:
        plotUpperLimits(labels, values)


if __name__ == "__main__":
    main()
