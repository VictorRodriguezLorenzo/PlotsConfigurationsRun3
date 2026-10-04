"""Unrolled 2D categorical-DNN discriminants for all Run-3 campaigns.

The inference aliases are expected to expose, for every mediator mass, three
consecutive softmax outputs in the order ``background, ttDM, tWDM``.  ROOT's
``y:x`` convention is used below, so the horizontal axis is the ttDM score and
the vertical axis is the tWDM score.  mkShapesRDF stores the TH2 and unfolds it
to a one-dimensional template when producing the datacard.
"""

blindSR = True
blindCuts = {
    f"ttdm_sr_{category}": "full"
    for category in cuts["ttdm_sr"]["categories"]  # noqa: F821 - injected by mkShapesRDF
} if blindSR else {}

variables = {}

# Coarse low-score bins retain background control, while the finer high-score
# bins preserve separation in the signal-rich corners without making 100-bin
# templates.  The same binning is used in every campaign for straightforward
# Run-3 merging.
scoreBins = [0.0, 0.05, 0.15, 0.30, 0.50, 0.70, 0.85, 0.95, 1.0]
mPhi = [50, 100, 150, 200, 250, 300, 350, 400, 500, 600, 700, 800, 1000]
classIndex = {"background": 0, "ttDM": 1, "tWDM": 2}

for category in ("ttDM", "tWDM"):
    for model in ("ps", "s"):
        alias = f"evaluate_dnn_categorical_{category}_{model}"
        for massIndex, mass in enumerate(mPhi):
            ttIndex = 3 * massIndex + classIndex["ttDM"]
            twIndex = 3 * massIndex + classIndex["tWDM"]
            variables[f"categorical_2D_{category}_{model}_{mass}"] = {
                "name": f"{alias}[{twIndex}]:{alias}[{ttIndex}]",
                "range": (scoreBins, scoreBins),
                "xaxis": (
                    f"categorical DNN p(ttDM):p(tWDM), {category}, "
                    f"{model}, m_{{#Phi}} = {mass} GeV"
                ),
                "fold": 3,
                "blind": blindCuts,
            }
