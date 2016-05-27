clc
clear all
close all

n = linspace(1, 6, 100);
LIMIT_FOS = 1;

FoS_incTorsion = zeros(size(n));
FoSUlt_incTorsion = zeros(size(n));
FoSSkin = zeros(size(n));
FoSSkin_curv = zeros(size(n));
FoSir = zeros(size(n));
FoScr = zeros(size(n));
FoScripPure = zeros(size(n));
FoScripAdjusted = zeros(size(n));
FoSStringer = zeros(size(n));

LIMIT_LOAD = 3;
ULT_LOAD = 5;

for i = 1:length(n)
    [ FoS_incTorsion(i), FoSUlt_incTorsion(i), FoSSkin(i), FoSSkin_curv(i), FoSir(i), FoSCr(i), FoScripPure(i), FoScripAdjusted(i), FoSStringer(i) ] ...
        = analysisTotal(n(i), 0);
end

figure;
plot(n, FoSir, n, FoS_incTorsion, n, FoSSkin_curv, n, FoSCr, n, FoScripAdjusted, n, ones(size(n)), 'k--');
xlabel('Load Factor');
ylabel('Factor of Safety');
legend('Interfastener Buckling', 'Yielding', 'Skin Panel Buckling (w/ curvature)', 'Stiffener Column Buckling', 'Crippling' ,'Failure Margin')
title('Failure modes of wing section with load factor');

idxYield = find(FoS_incTorsion < LIMIT_FOS, 1, 'first');
LFYield = n(idxYield)
idxSkin = find(FoSSkin < LIMIT_FOS, 1, 'first');
LFSkin = n(idxSkin)
idxCr = find(FoSCr < LIMIT_FOS, 1, 'first');
LFCr = n(idxCr)
idxCrip = find(FoScripAdjusted < LIMIT_FOS, 1, 'first');
LFCrip = n(idxCrip)

figure;
plot(n, FoSSkin, n, FoSSkin_curv, n, ones(size(n)), 'k--');
xlabel('Load Factor');
ylabel('Factor of Safety');
legend('Skin Panel Buckling (no curvature)', 'Skin Panel Buckling (w/ curvature)');
title('Effect of curvature on skin panel buckling')

