import ROOT
from ROOT import TCanvas, TGraph, TLegend

DEFAULT_LUMI_LABEL = "109 fb^{-1} (13.6 TeV)"
ADD_OBSERVED = True


def getLimits(filename):
    limits = []
    with open(filename) as f:
        for line in f:
            if "Observed Limit" in line:
                limits.append(float(line.split()[-1]))
            elif "Expected  2.5%" in line:
                limits.append(float(line.split()[-1]))
            elif "Expected 16.0%" in line:
                limits.append(float(line.split()[-1]))
            elif "Expected 50.0%" in line:
                limits.append(float(line.split()[-1]))
            elif "Expected 84.0%" in line:
                limits.append(float(line.split()[-1]))
            elif "Expected 97.5%" in line:
                limits.append(float(line.split()[-1]))
    return limits


def plotUpperLimits(labels, values, lumi_label=DEFAULT_LUMI_LABEL, model="ps", output_dir="."):
    lumi_label = lumi_label or DEFAULT_LUMI_LABEL
    ROOT.gROOT.SetBatch(True)

    if model not in {"ps", "s"}:
        raise ValueError(f"Unsupported model {model!r}; expected 'ps' or 's'")

    N = len(labels)
    yellow = TGraph(2 * N)
    green = TGraph(2 * N)
    median = TGraph(N)
    observed = TGraph(N)

    up2s = []
    draw_observed = ADD_OBSERVED

    for i in range(N):
        file_name = labels[i]
        limit = getLimits(file_name)

        if len(limit) == 6:
            observed_limit = limit[0]
            expected = limit[1:]
        elif len(limit) == 5:
            observed_limit = None
            expected = limit
            draw_observed = False
        else:
            raise ValueError(
                f"Expected 5 expected-limit values, optionally with 1 observed value, in {file_name}; found {len(limit)}"
            )

        up2s.append(expected[4])

        yellow.SetPoint(i, values[i], expected[4])
        green.SetPoint(i, values[i], expected[3])
        median.SetPoint(i, values[i], expected[2])
        green.SetPoint(2 * N - 1 - i, values[i], expected[1])
        yellow.SetPoint(2 * N - 1 - i, values[i], expected[0])

        if draw_observed and observed_limit is not None:
            observed.SetPoint(i, values[i], observed_limit)

    W = 750
    H = 800
    T = 0.08 * H
    B = 0.12 * H
    L = 0.12 * W
    R = 0.04 * W

    c = TCanvas(f"c_{model}", f"c_{model}", 100, 100, W, H)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)
    c.SetLeftMargin(L / W)
    c.SetRightMargin(R / W)
    c.SetTopMargin(T / H)
    c.SetBottomMargin(B / H)
    c.SetTickx(0)
    c.SetTicky(0)
    c.SetGrid(0, 0)
    c.SetLogy()
    c.cd()

    x_min = min(values) * 0.9
    x_max = max(values) * 1.1
    y_min = 1e-2
    y_max = max(up2s) * 1.6

    frame = c.DrawFrame(x_min, y_min, x_max, y_max)
    frame.GetYaxis().CenterTitle()
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetXaxis().SetTitleOffset(1.0)
    frame.GetYaxis().SetTitleOffset(0.9)
    frame.GetXaxis().SetNdivisions(508)
    frame.GetYaxis().CenterTitle(True)
    frame.GetYaxis().SetTitle("95% CL upper limit on #mu = #sigma / #sigma_{theory}")
    frame.GetXaxis().SetTitle("m_{a} [GeV]" if model == "ps" else "m_{#phi} [GeV]")
    frame.SetMinimum(y_min)
    frame.SetMaximum(y_max)
    frame.GetXaxis().SetLimits(10, 1000)

    green_col = ROOT.TColor.GetColor("#6B8E4E")
    yellow_col = ROOT.TColor.GetColor("#F6C35B")

    yellow.SetFillColor(yellow_col)
    yellow.SetLineColor(yellow_col)
    yellow.SetLineWidth(0)
    yellow.SetFillStyle(1001)
    yellow.Draw("F")

    green.SetFillColor(green_col)
    green.SetLineColor(green_col)
    green.SetLineWidth(0)
    green.SetFillStyle(1001)
    green.Draw("Fsame")

    median.SetLineColor(ROOT.kBlack)
    median.SetLineWidth(2)
    median.SetLineStyle(2)
    median.Draw("Lsame")

    if draw_observed:
        observed.SetLineColor(ROOT.kBlack)
        observed.SetLineWidth(2)
        observed.SetMarkerStyle(20)
        observed.SetMarkerSize(1.0)
        observed.Draw("LPsame")

    ref_line = ROOT.TLine(frame.GetXaxis().GetXmin(), 1.0, frame.GetXaxis().GetXmax(), 1.0)
    ref_line.SetLineColor(ROOT.kRed + 1)
    ref_line.SetLineStyle(2)
    ref_line.SetLineWidth(2)
    ref_line.Draw("same")

    cms = ROOT.TLatex()
    cms.SetNDC()
    cms.SetTextFont(61)
    cms.SetTextSize(0.06)
    cms.DrawLatex(0.13, 0.93, "CMS")

    pre = ROOT.TLatex()
    pre.SetNDC()
    pre.SetTextFont(52)
    pre.SetTextSize(0.045)
    pre.DrawLatex(0.28, 0.93, "Preliminary")

    lumi = ROOT.TLatex()
    lumi.SetNDC()
    lumi.SetTextFont(42)
    lumi.SetTextSize(0.045)
    lumi.SetTextAlign(31)
    lumi.DrawLatex(0.95, 0.93, lumi_label)

    model_text = ROOT.TLatex()
    model_text.SetNDC()
    model_text.SetTextFont(42)
    model_text.SetTextAlign(11)
    model_text.SetTextSize(0.040)
    if model == "ps":
        model_text.DrawLatex(0.52, 0.37, "Pseudoscalar a #rightarrow #chi#bar{#chi}")
    else:
        model_text.DrawLatex(0.52, 0.37, "Scalar #phi #rightarrow #chi#bar{#chi}")
    model_text.SetTextSize(0.036)
    model_text.DrawLatex(0.52, 0.32, "g_{q} = g_{#chi} = 1")
    model_text.DrawLatex(0.52, 0.27, "Dirac DM, m_{#chi} = 1 GeV")

    leg_green = TGraph()
    leg_green.SetFillColor(green_col)
    leg_green.SetLineColor(green_col)
    leg_green.SetLineWidth(0)
    leg_green.SetFillStyle(1001)

    leg_yellow = TGraph()
    leg_yellow.SetFillColor(yellow_col)
    leg_yellow.SetLineColor(yellow_col)
    leg_yellow.SetLineWidth(0)
    leg_yellow.SetFillStyle(1001)

    legend = TLegend(0.19, 0.67, 0.47, 0.85)
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextFont(42)
    legend.SetTextSize(0.038)
    legend.AddEntry(median, "Expected", "L")
    legend.AddEntry(leg_green, "#pm1 #sigma expected", "F")
    legend.AddEntry(leg_yellow, "#pm2 #sigma expected", "F")
    if draw_observed:
        legend.AddEntry(observed, "Observed", "LP")
    legend.Draw()

    ROOT.gPad.SetTicks(1, 1)
    frame.Draw("sameaxis")

    import os
    os.makedirs(output_dir, exist_ok=True)
    stem = f"BrazilBplot_{model}"
    c.SaveAs(os.path.join(output_dir, stem + ".png"))
    c.SaveAs(os.path.join(output_dir, stem + ".pdf"))
    c.Close()
