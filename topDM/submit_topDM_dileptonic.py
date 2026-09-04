import subprocess
import glob
import os

# Mediator masses
mPhi = ['10','50','100','150','200','250','300','350','400','500','600','700','800','1000']

samples = []
for phi in mPhi:
    samples.append(f"TTto2LDMsimpSpin0_s_mphi-{phi}")
    samples.append(f"TTto2LDMsimpSpin0_ps_mphi-{phi}")
    samples.append(f"TTto2LDMsimpSpin0_s_mphi-{phi}_ext1")
    samples.append(f"TTto2LDMsimpSpin0_ps_mphi-{phi}_ext1")
    samples.append(f"TWto2LDMsimpSpin0_s_mphi-{phi}")
    samples.append(f"TWto2LDMsimpSpin0_ps_mphi-{phi}")

campaign_bases = {
#    "Summer22_130x_nAODv12_Full2022v12": "MCl2loose2022v12__MCCorr2022v12JetScaling__l2tight",
#    "Summer22EE_130x_nAODv12_Full2022v12": "MCl2loose2022EEv12__MCCorr2022EEv12JetScaling__l2tight",
#    "Summer23_130x_nAODv12_Full2023v12": "MCl2loose2023v12__MCCorr2023v12JetScaling__l2tight",
#    "Summer23BPix_130x_nAODv12_Full2023BPixv12": "MCl2loose2023BPixv12__MCCorr2023BPixv12JetScaling__l2tight",
    "Summer24_150x_nAODv15_Full2024v15": "MCl2loose2024v15__MCCorr2024v15__JERFrom23BPix__l2tight"
}

for campaign, base in campaign_bases.items():
    base_dir = f"/eos/user/v/victorr/HWWNano/{campaign}/{base}"

    for s in samples:
        pattern = os.path.join(base_dir, f"nanoLatino_{s}__part*.root")

        cmd = [
            "mkPostProc",
            "-o", "0",
            "-p", campaign,
            "-s", base,
            "-T", s
        ]

        print("\nRunning:", " ".join(cmd))
        result = subprocess.run(cmd)

        # Retry with --useRedirector if it crashed
        if result.returncode != 0:
            print(f"ERROR running {s} in {campaign}, retrying with --useRedirector...\n")
            cmd.append("--useRedirector 1")
            result_retry = subprocess.run(cmd)
            if result_retry.returncode != 0:
                print(f"FAILED again for {s} in {campaign}, skipping...\n")
            else:
                print(f"Success on retry for {s} with --useRedirector.\n")
