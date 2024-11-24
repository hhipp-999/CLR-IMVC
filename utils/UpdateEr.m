function [Er] = UpdateEr(H,Zr,Er,Y4,mu,lambda4)
        F = [];
        F = [F,H-H*Zr+Y4/mu];
        [Econcat] = solve_l1l2(F,lambda4/mu);
        Er = [Econcat];