% Function that returns the value of the B5D36 spectrometer dipole value in
% GeV, if it was saved and if it is not a nan. Otherwise, it returns the
% nominal 20.35 GeV setting.

function B5D36 = getB5D36(E200_state)

if isfield(E200_state, 'LI20_LGPS_3330_BDES')
    if ~isnan(E200_state.LI20_LGPS_3330_BDES);
        B5D36 = E200_state.LI20_LGPS_3330_BDES; 
    else
        B5D36 = 20.35;
        disp('Not possible to access B5D36 value in E200_state. B5D36 set to nominal 20.35 GeV.');
    end
else
    B5D36 = 20.35;
    disp('E200_state was not saved. B5D36 set to nominal 20.35 GeV.');
end
end