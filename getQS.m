% Function that returns the value of the QS energy setting in
% GeV (relative to 20.35 GeV), if it was saved, if QS1 and QS2 are not nan
% and if the values correspond to a QS energy setting that images plasma
% exit in 2013. Otherwise, it returns the nominal 0 GeV setting.


function QS_val = getQS(E200_state)

if isfield(E200_state, 'LI20_LGPS_3261_BDES')
    if ~isnan(E200_state.LI20_LGPS_3261_BDES) && ~isnan(E200_state.LI20_LGPS_3311_BDES);
        QS = [E200_state.LI20_LGPS_3261_BDES, E200_state.LI20_LGPS_3311_BDES];
        E = 20.35 * QS./[213.07, -156.01] - 20.35;
        if abs(E(1)-E(2))<0.1
            QS_val = E(1);
        else
            QS_val = 0;
            disp('QS 1 and QS 2 values does not match a QS energy setting that images plasma exit in 2013. QS_val set to nominal 0 GeV.')
        end
    else
        QS_val = 0;
        disp('Not possible to access QS1 or QS2 value in E200_state. QS_val set to nominal 0 GeV.')
    end
else
    QS_val = 0;
    disp('E200_state was not saved. QS_val set to nominal 0 GeV.');
end
end