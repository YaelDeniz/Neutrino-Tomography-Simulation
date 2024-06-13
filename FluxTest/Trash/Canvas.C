void Canvas(){
   const double x[] = {0.0, 1.0};  
   const double y[] = {20000.0, 10000.0};
   TGraph g(2, x, y);
   g.GetYaxis() -> SetTitle("#frac{A}{B}");
   g.GetYaxis() -> SetTitleOffset(1.);
   g.Draw("A*");
}