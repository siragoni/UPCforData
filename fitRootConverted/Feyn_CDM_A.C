// draw feynman diagram of exclusive photo production

void Feyn_CDM_A()
{

   TCanvas *c2 = new TCanvas("c2", "c2", 10,10, 600, 400);
   c2->Range(0, 0, 100, 80);
   Int_t linsav = gStyle->GetLineWidth();
   gStyle->SetLineWidth(2);
   TLatex t;
   t.SetTextAlign(22);
   t.SetTextSize(0.05);
   TLine * l;

   //draw arc
   TArc *a = new TArc(50,20,5,195,344); a->SetLineStyle(2);
   a->Draw();
   t.DrawLatex(50,10,"t");

   TArc *a1 = new TArc(58,29,25,150,210); a1->SetLineStyle(3);
   a1->SetNoEdges();
   a1->Draw();
   //  t.DrawLatex(28,30,"W_{#gammap}");
     t.DrawLatex(28,30,"W");


   // draw upper lead
   l = new TLine(10, 70, 25, 60); l->Draw();
   l = new TLine(25, 60, 90, 70); l->Draw();
   //   t.DrawLatex(15,74,"p, Pb");
   // t.DrawLatex(85,74,"p, Pb");
   t.DrawLatex(15,74,"A");
   t.DrawLatex(85,74,"A");

   // draw lower particle
   l = new TLine(10, 10, 50, 20); l->Draw();
   l = new TLine(50, 20, 90, 10); l->Draw();
   t.DrawLatex(15,6,"p, A");
   t.DrawLatex(85,6,"p, A");

   // draw photon
   TCurlyLine *gamma = new TCurlyLine(25, 60, 40, 37);
   gamma->SetWavy();
   gamma->Draw();
   t.DrawLatex(28,47,"#gamma");

   //draw dipole
   TEllipse *c = new TEllipse(50,37,12,4,90,-90,180);
   c->SetNoEdges();
   c->Draw();
   TEllipse *c1 = new TEllipse(50,37,12,4,90,-90);
   c1->SetNoEdges();
   c1->Draw();


   // draw VM
   l = new TLine(62, 37, 85, 37); l->Draw();
   t.DrawLatex(73,47,"Vector meson");
   t.DrawLatex(73,42,"(#rho^{0}, J/#psi, #psi(2S), ...)");
   //  t.DrawLatex(73,44,"J/#psi");

   // draw interaction
   TEllipse *i = new TEllipse(50,30,14,3,0,360,90);
   i->SetNoEdges();
   i->SetLineColor(17);
   i->SetFillColor(17);
   i->Draw();



}


//_____________________________________________________________________________
/* - Feynman diagram for the J/Psi.
 * -
 */
void Feyn_CDM_B()
{

   TCanvas *c2 = new TCanvas("c2", "c2", 10,10, 600, 400);
   c2->Range(0, 0, 100, 80);
   Int_t linsav = gStyle->GetLineWidth();
   gStyle->SetLineWidth(2);
   TLatex t;
   t.SetTextAlign(22);
   t.SetTextSize(0.05);
   t.SetTextSize(0.07);

   TLine * l;
   gStyle->SetLineWidth(2);

   //draw arc
   TArc *a = new TArc(50,20,5,195,344); a->SetLineStyle(2);
   a->Draw();
   t.DrawLatex(50,12,"t");
   // t.DrawLatex(50, 3,"(a)");

   TArc *a1 = new TArc(58,29,25,150,210); a1->SetLineStyle(3);
   a1->SetNoEdges();
   a1->Draw();
   //  t.DrawLatex(28,30,"W_{#gammap}");
     t.DrawLatex(28,30,"W");


   // draw upper lead
   l = new TLine(10, 70, 25, 60); l->Draw();
   l = new TLine(25, 60, 90, 70); l->Draw();
   //   t.DrawLatex(15,74,"p, Pb");
   // t.DrawLatex(85,74,"p, Pb");
   t.DrawLatex(15,74,"Pb");
   t.DrawLatex(85,74,"Pb");

   // draw lower particle
   l = new TLine(10, 10, 50, 20); l->Draw();
   l = new TLine(50, 20, 90, 10); l->Draw();
   t.DrawLatex(15,6,"p, Pb");
   t.DrawLatex(85,6,"p, Pb");

   // draw photon
   TCurlyLine *gamma = new TCurlyLine(25, 60, 40, 37);
   gamma->SetWavy();
   gamma->Draw();
   t.DrawLatex(28,47,"#gamma");

   //draw dipole
   TEllipse *c = new TEllipse(50,37,12,4,90,-90,180);
   c->SetNoEdges();
   c->Draw();
   TEllipse *c1 = new TEllipse(50,37,12,4,90,-90);
   c1->SetNoEdges();
   c1->Draw();


   // draw VM
   l = new TLine(62, 37, 85, 37); l->Draw();
   t.SetTextSize(0.07);
   t.DrawLatex(76,47,"Vector meson");
   t.DrawLatex(76,42,"(J/#psi, #psi(2S), ...)");
   //  t.DrawLatex(73,44,"J/#psi");

   // draw interaction
   TEllipse *i = new TEllipse(50,30,14,3,0,360,90);
   i->SetNoEdges();
   i->SetLineColor(17);
   i->SetFillColor(17);
   i->Draw();



}

//_____________________________________________________________________________
/* - Feynman diagram for the J/Psi.
 * -
 */
void Feyn_CDM_B2()
{

   TCanvas *c2 = new TCanvas("c2", "c2", 10,10, 600, 400);
   c2->Range(0, 0, 100, 80);
   Int_t linsav = gStyle->GetLineWidth();
   gStyle->SetLineWidth(2);
   TLatex t;
   t.SetTextAlign(22);
   t.SetTextSize(0.05);
   t.SetTextSize(0.07);

   TLine * l;
   gStyle->SetLineWidth(2);

   //draw arc
   // TArc *a = new TArc(50,20,5,195,344); a->SetLineStyle(2);
   // a->Draw();
   // t.DrawLatex(50,12,"t");
   // t.DrawLatex(50, 3,"(b)");

   // TArc *a1 = new TArc(70,29,25,150,210); a1->SetLineStyle(3);
   // a1->SetNoEdges();
   // a1->Draw();
   // //  t.DrawLatex(28,30,"W_{#gammap}");
   //   t.DrawLatex(28,47,"W");


   // draw upper lead
   l = new TLine(10, 70, 50, 60); l->Draw();
   l = new TLine(50, 60, 90, 70); l->Draw();
   //   t.DrawLatex(15,74,"p, Pb");
   // t.DrawLatex(85,74,"p, Pb");
   t.DrawLatex(15,74,"Pb, p");
   t.DrawLatex(85,74,"Pb, p");

   // draw lower particle
   l = new TLine(10, 10, 25, 20); l->Draw();
   l = new TLine(25, 20, 90, 10); l->Draw();
   t.DrawLatex(15,6,"Pb, p");
   t.DrawLatex(85,6,"Pb, p");

   // draw photon
   TCurlyLine *gamma = new TCurlyLine(25, 20, 40, 35);
   gamma->SetWavy();
   gamma->Draw();
   t.DrawLatex(28,30,"#gamma");

   //draw dipole
   TEllipse *c = new TEllipse(50,37,12,4,90,-90,180);
   c->SetNoEdges();
   c->Draw();
   TEllipse *c1 = new TEllipse(50,37,12,4,90,-90);
   c1->SetNoEdges();
   c1->Draw();


   // draw VM
   l = new TLine(62, 37, 85, 37); l->Draw();
   t.SetTextSize(0.07);
   t.DrawLatex(76,47,"Vector meson");
   t.DrawLatex(76,42,"(J/#psi, #psi(2S), ...)");
   //  t.DrawLatex(73,44,"J/#psi");

   // draw interaction
   TEllipse *i = new TEllipse(50,45,17,5,0,360,90);
   i->SetNoEdges();
   i->SetLineColor(17);
   i->SetFillColor(17);
   i->Draw();



}
//_____________________________________________________________________________
/* - Feynman diagram for the GammaGammaBkg.
 * -
 */
void Feyn_CDM_C()
{

   TCanvas *c2 = new TCanvas("c2", "c2", 10,10, 600, 400);
   c2->Range(0, 0, 100, 80);
   Int_t linsav = gStyle->GetLineWidth();
   gStyle->SetLineWidth(2);
   TLatex t;
   t.SetTextAlign(22);
   t.SetTextSize(0.05);
   t.SetTextSize(0.07);

   TLine * l;
   gStyle->SetLineWidth(2);

   t.DrawLatex(50, 3,"(b)");



   // draw upper lead
   l = new TLine(10, 70, 25, 60); l->Draw();
   l = new TLine(25, 60, 90, 70); l->Draw();
   //   t.DrawLatex(15,74,"p, Pb");
   // t.DrawLatex(85,74,"p, Pb");
   t.DrawLatex(15,74,"Pb, p");
   t.DrawLatex(85,74,"Pb, p");

   // draw lower particle
   l = new TLine(10, 10, 25, 20); l->Draw();
   l = new TLine(25, 20, 90, 10); l->Draw();
   t.DrawLatex(15,6,"Pb, p");
   t.DrawLatex(85,6,"Pb, p");

   // draw photon
   TCurlyLine *gamma  = new TCurlyLine(25, 60, 45, 45);
   gamma->SetWavy();
   gamma->Draw();
   // t.DrawLatex(28,47,"#gamma");

   // draw second photon
   TCurlyLine *gamma2 = new TCurlyLine(25, 20, 45, 35);
   gamma2->SetWavy();
   gamma2->Draw();
   t.DrawLatex(28,50,"#gamma");
   t.DrawLatex(28,30,"#gamma");


   // draw virtual lepton
   l = new TLine(45, 35, 45, 45); l->Draw();
   l = new TLine(45, 35, 85, 30); l->Draw();
   l = new TLine(45, 45, 85, 50); l->Draw();
   t.DrawLatex(76,53,"#mu^{+}");
   t.DrawLatex(76,27,"#mu^{-}");




}
//_____________________________________________________________________________
//_____________________________________________________________________________
/* - Feynman diagram for the J/Psi.
 * -
 */
void Feyn_CDM_pPb()
{

   TCanvas *c2 = new TCanvas("c2", "c2", 10,10, 600, 400);
   c2->Range(0, 0, 100, 80);
   Int_t linsav = gStyle->GetLineWidth();
   gStyle->SetLineWidth(2);
   TLatex t;
   t.SetTextAlign(22);
   t.SetTextSize(0.05);
   t.SetTextSize(0.07);

   TLine * l;
   gStyle->SetLineWidth(2);

   //draw arc
   TArc *a = new TArc(50,20,5,195,344); a->SetLineStyle(2);
   a->Draw();
   t.DrawLatex(50,12,"t");
   // t.DrawLatex(50, 3,"(a)");

   TArc *a1 = new TArc(58,29,25,150,210); a1->SetLineStyle(3);
   a1->SetNoEdges();
   a1->Draw();
   //  t.DrawLatex(28,30,"W_{#gammap}");
     t.DrawLatex(28,30,"W");


   // draw upper lead
   l = new TLine(10, 70, 25, 60); l->Draw();
   l = new TLine(25, 60, 90, 70); l->Draw();
   //   t.DrawLatex(15,74,"p, Pb");
   // t.DrawLatex(85,74,"p, Pb");
   t.DrawLatex(15,74,"Pb");
   t.DrawLatex(85,74,"Pb");

   // draw lower particle
   l = new TLine(10, 10, 50, 20); l->Draw();
   l = new TLine(50, 20, 90, 10); l->Draw();
   t.DrawLatex(15,6,"p");
   t.DrawLatex(85,6,"p");

   // draw photon
   TCurlyLine *gamma = new TCurlyLine(25, 60, 40, 37);
   gamma->SetWavy();
   gamma->Draw();
   t.DrawLatex(28,47,"#gamma");

   //draw dipole
   TEllipse *c = new TEllipse(50,37,12,4,90,-90,180);
   c->SetNoEdges();
   c->Draw();
   TEllipse *c1 = new TEllipse(50,37,12,4,90,-90);
   c1->SetNoEdges();
   c1->Draw();


   // draw VM
   l = new TLine(62, 37, 85, 37); l->Draw();
   t.SetTextSize(0.07);
   t.DrawLatex(76,47,"Vector meson");
   t.DrawLatex(76,42,"(J/#psi, #psi(2S), ...)");
   //  t.DrawLatex(73,44,"J/#psi");

   // draw interaction
   TEllipse *i = new TEllipse(50,30,14,3,0,360,90);
   i->SetNoEdges();
   i->SetLineColor(17);
   i->SetFillColor(17);
   i->Draw();



}
//_____________________________________________________________________________
/* - Feynman diagram for the GammaGammaBkg.
 * -
 */
void FullDraw()
{
  TCanvas* Can = new TCanvas("Can", "Can", 1200, 800);
  Can->Divide(2,1);
  Can->cd(1);
  Feyn_CDM_B();
  Can->cd(2);
  Feyn_CDM_C();
}
