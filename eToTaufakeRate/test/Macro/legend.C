

{




  TCanvas *c1 = new TCanvas("c1","Canvas",10,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TH1F* h = new TH1F("h","",2,0,1);

  h->SetMarkerColor(kBlue);
  h->SetMarkerStyle(kOpenCircle);
  h->SetMarkerSize(1.5);

  TLegend* leg = new TLegend(0.4,0.15,0.85,0.5,NULL,"brNDC");
  leg->SetFillStyle(4000);
  leg->SetBorderSize(0);
  //leg->SetFillColor(1);
  leg->SetTextSize(0.03);

  leg->AddEntry(h,"ciao mondo","LEP");

  h->Draw();
  leg->Draw();

}
