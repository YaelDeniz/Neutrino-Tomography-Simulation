%% Energy Resolution
EOmin = 1.5
EOmax = 10
Ae=0.2
Be=0.0
myfun = @(Emin) -Emin + EOmin - 4*(Ae*Emin + Be*sqrt(Emin))
myfun2 = @(Emax)  -Emax + EOmax + 4*(Ae*Emax + Be*sqrt(Emax))
E0 = 1;
Elow = fzero(myfun,1) %Minimum Energy
Ehigh = fzero(myfun2,190)

syms x
fplot(myfun2(x), [0 50])

%% Angular Resolution
thOmin = 140
thOmax = 170
Ath=2.0
Bth=10.0
myfun = @(thmin) -thmin + thOmin - 4*(Ath + Bth/sqrt(Elow))
myfun2 = @(thmax)  -thmax + thOmax + 4*(Ath + Bth/sqrt(Elow))
E0 = 180;
E1 = fzero(myfun,1)
E2 = fzero(myfun2,190)

syms x
fplot(myfun2(x), [0 50])
