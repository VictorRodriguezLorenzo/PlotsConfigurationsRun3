import os

user = os.environ["USER"]

foldersToMerge = {
    "Full2022v12": {
        "folder": "../Full2022v12",
        "tag": "ttDM_dilep_2022",
        "outputFolder": f"/eos/user/{user[0]}/{user}/www/Run3-ttDM/rootFiles/ttDM_dilep_2022/",
    },
    "Full2022EEv12": {
        "folder": "../Full2022EEv12",
        "tag": "ttDM_dilep_2022EE",
        "outputFolder": f"/eos/user/{user[0]}/{user}/www/Run3-ttDM/rootFiles/ttDM_dilep_2022EE/",
    },
    "Full2023v12": {
        "folder": "../Full2023v12",
        "tag": "ttDM_dilep_2023",
        "outputFolder": f"/eos/user/{user[0]}/{user}/www/Run3-ttDM/rootFiles/ttDM_dilep_2023/",
    },
    "Full2023BPixv12": {
        "folder": "../Full2023BPixv12",
        "tag": "ttDM_dilep_2023BPix",
        "outputFolder": f"/eos/user/{user[0]}/{user}/www/Run3-ttDM/rootFiles/ttDM_dilep_2023BPix/",
    },
    "Full2024v15": {
        "folder": "../Full2024v15",
        "tag": "ttDM_dilep_2024",
        "outputFolder": f"/eos/user/{user[0]}/{user}/www/Run3-ttDM/rootFiles/ttDM_dilep_2024/",
    },
}


#    "Full2025v15": {
#        "folder": "../Full2025v15",
#        "tag": "ttDM_dilep_2025",
#        "outputFolder": f"/eos/user/{user[0]}/{user}/www/Run3-ttDM/rootFiles/ttDM_dilep_2025/",
#    },
