%a function that knowing t, ts, As/bt, and ts/t, calculate bs
%bs = wb.t_Upper.*(sigmaY./(wb.sr_s0.*Ks.*E)).^(-1/2)
function [bs,hs,Ks,F,sr_s0] = findbs(t, sigmaY, srs0_guess, Ks_guess, E,As_bt,ts_t,tol)
    Ks = Ks_guess.*ones(size(t));
    sr_s0 = srs0_guess;
    ts_t = ts_t.*ones(size(ts_t));
    bs2 = 100.*ones(size(t));
    converged = false;
    maxiter = 100;
    iter = 1;
    while ~converged && iter < maxiter
        bs = max(0, t.*(sigmaY./(sr_s0.*Ks.*E)).^(-1/2)); %Stringer pitch
        %find hs from As/bt
        hs = As_bt./(ts_t.*1.6./(bs));
        Ks = interpKsSkinStringer(hs./bs, ts_t);
        As_bt = 1.6*hs./bs.*ts_t;
        sr_s0 = interpsr_s0(As_bt, ts_t);
        F = interpF(As_bt, ts_t);
        err = sqrt(nanmean((bs - bs2).^2));
        if err < tol
            converged = true;
        end
        bs2 = bs;
        iter = iter + 1;
    end
end