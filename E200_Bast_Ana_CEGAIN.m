% Script to analyze butterfly from E200 2013 data with butterfly on CEGAIN camera

% Bastien CLAUSSE

% Create: June 16, 2013
% Last edit: July 25, 2013

%% PATH SETTINGS

addpath('C:/Users/Bastien/Desktop/PRE/Matlab/E200_scripts/facet_daq');
addpath('C:/Users/Bastien/Desktop/PRE/Matlab/E200_scripts/bast_ana');

prefix = 'C:/Users/Bastien/Desktop/PRE/Matlab/Volumes/PWFA_4big';
day = '20130429';
data_set = 'E200_10854';
do_save = 0;
save_path = ['C:/Users/Bastien/Desktop/PRE/Matlab/2013_E200_Data_Analysis/Rubidium_density=6.6e16_fullpyro/Zoom_Butterfly_setQSenergy=-8GeV/' day '/'];


%% SCAN INFO FILE INITIALISATION

cmap = custom_cmap();

CEGAIN_caxis = [0.8 3.2];

path = [prefix '/nas/nas-li20-pm01/E200/2013/' day '/' data_set '/'];
if(~exist([save_path data_set], 'dir')); mkdir([save_path data_set '/frames/']); end;

scan_info_file = dir([path '*scan_info*']);
if size(scan_info_file,1) == 1
    load([path scan_info_file.name]);
    n_step = size(scan_info,2);
    is_qsbend_scan = strcmp(scan_info(1).Control_PV_name, 'set_QSBEND_energy');
    is_qs_scan = strcmp(scan_info(1).Control_PV_name, 'set_QS_energy');
    is_scan = 1;
elseif size(scan_info_file,1) == 0
    filenames_file = dir([path data_set '*_filenames.mat']);
    load([path filenames_file.name]);
    scan_info = filenames;
    n_step = 1;
    is_qsbend_scan = 0;
    is_qs_scan = 0;
    is_scan = 0;
else
    error('There are more than 1 scan info file.');
end

list = dir([path data_set '_2013*.mat']);
mat_filenames = {list.name};
mat_filenames = {mat_filenames{1:2:end}};
load([path mat_filenames{1}]);
n_shot = param.n_shot;

CEGAIN = cam_back.CEGAIN;

CEGAIN.X_RTCL_CTR = 700;
CEGAIN.Y_RTCL_CTR = 500;

CEGAIN.xx = 1e-3*CEGAIN.RESOLUTION * ( (CEGAIN.ROI_X-CEGAIN.X_RTCL_CTR+1):(CEGAIN.ROI_X+CEGAIN.ROI_XNP-CEGAIN.X_RTCL_CTR) );
CEGAIN.yy = 1e-3*CEGAIN.RESOLUTION * ( (CEGAIN.ROI_Y-CEGAIN.Y_RTCL_CTR+1):(CEGAIN.ROI_Y+CEGAIN.ROI_YNP-CEGAIN.Y_RTCL_CTR) );

clear processed;

processed.scalars.DIVERGENCE = zeros(n_step, n_shot);
processed.scalars.SIGMA0 = zeros(n_step, n_shot);
processed.scalars.y0 = zeros(n_step, n_shot);
processed.scalars.EMITTANCE = zeros(n_step, n_shot);

B5D36 = getB5D36(E200_state);
QS = getQS(E200_state);

E_EGAIN = E200_cher_get_E_axis('20130423', 'CEGAIN', 0, 1:1392, 0, B5D36);

order_fit = 2;

%% FIGURE SETTINGS

fig = figure(1);       
set(fig, 'position', [20, 50, 1340, 630]);
set(fig, 'PaperPosition', [0.25, 2.5, 35, 17]);
set(fig, 'color', 'w');
clf();

%% REQUIRED CALCULATIONS FOR STUDY

% for i=1:n_step    
i=4;

data = load([path mat_filenames{i}]);
if is_qsbend_scan; B5D36 = 20.35 + scan_info(i).Control_PV; end;

[CEGAIN.img, ~, CEGAIN.pid] = E200_readImages([prefix scan_info(i).CEGAIN]);

CEGAIN.img = double(CEGAIN.img);

for j=1:size(CEGAIN.img,3); CEGAIN.img(:,:,j) = CEGAIN.img(:,:,j) - cam_back.CEGAIN.img(:,:); end;

%% ZOOM SETTINGS

    nb_slice = 30; % it represents the number of slices, and so the number of points you will use for the fit
                   % keep in mind it's better to conserv (xzoom_max - xzoom_min + 1)/nb_slice (and so the number of pixels per slice) integer
   
    [~,ind]=ismember(min((E_EGAIN-(B5D36 + scan_info(i).Control_PV)).^2),(E_EGAIN-(B5D36 + scan_info(i).Control_PV)).^2);
    xzoom_min = ind -90;
    xzoom_max = ind +89;
 
    yzoom_min = 350;
    yzoom_max = 625;
 
    
%%  
  for j=1:n_shot
% j = 29;
%%
    fprintf('\n%d ', CEGAIN.pid(j));
    
    [CEGAIN.ana, CEGAIN.ana.img] = Ana_CEGAIN_img(E_EGAIN, CEGAIN.img(:,:,j));
            
    CEGAIN.ana.img(CEGAIN.ana.img<1) = 1;
    CEGAIN.img2 = CEGAIN.img(:,:,j);
    CEGAIN.img2(CEGAIN.img2<1) = 1;

%%  SLICES AND SIGMA  
  
    CEGAIN.yyzoom = CEGAIN.yy(yzoom_min:yzoom_max);
    CEGAIN.xxzoom = CEGAIN.xx(xzoom_min:xzoom_max);
    CEGAIN.zoom = CEGAIN.ana.img(yzoom_min:yzoom_max, xzoom_min:xzoom_max); % Using filters        
    CEGAIN.slice = slice(CEGAIN.zoom, nb_slice);
       
    h = floor((xzoom_max - xzoom_min +1)/nb_slice); % number of pixels per slice
    
    CEGAIN.reduc = reduc(CEGAIN.slice, nb_slice); % it defines the line-outs that will be plotted.
    CEGAIN.reduc1 = CEGAIN.reduc(:,h+1);
    CEGAIN.reduc2 = CEGAIN.reduc(:,h*0.5*nb_slice + 1);
    CEGAIN.reduc3 = CEGAIN.reduc(:,1+h*(nb_slice-1));
     
    yparabole2=[];
    for k=0:(nb_slice - 1); yparabole2 = [yparabole2, CEGAIN.xxzoom(1,1+k*h)];end;
    
    A2 = sigmacarre(CEGAIN.slice', CEGAIN.yyzoom, nb_slice); % A2 is a vector with the different values of sigma squared with the desired zoom.
           
    p2 = polyfit(yparabole2, A2,order_fit); % Fit calculation
    x3 = yparabole2;
    y3 = polyval(p2,x3);
    
    eta2 = 0.001*61*B5D36/(B5D36 + scan_info(i).Control_PV);   % in m
    theta2 = sqrt(0.25*p2(order_fit - 1)*(eta2^2)/(23.05^2));    % in rad
    y02 = -p2(order_fit)/(2*p2(order_fit - 1)); % in mm
    sigma02 =  sqrt(p2(order_fit + 1)-p2(order_fit - 1)*y02^2)/(5.6); % in mm
    emit2 = 40000*sigma02*1000*theta2; % in microns

%% PROCESSED
         
         processed.scalars.DIVERGENCE(i,j) = 1000000*theta2;  % in urad
         processed.scalars.SIGMA0(i,j) = 1000*sigma02;  % in um
         processed.scalars.y0(i,j) = y02;  % in mm
         processed.scalars.EMITTANCE(i,j) = emit2;  % in microns  
    
%%  LOOPS INITIALISATION 
    if i==3 && j==1  
%%   SETTINGS DISPLAY OF SHOT INFORMATION AND FITTED PARAMETERS 

   h_text = axes('Position', [0.07, 0.80, 0.2, 0.05], 'Visible', 'off');    
        if is_scan; STEP = text(0., 0., [char(regexprep(scan_info(i).Control_PV_name, '_', '\\_')) ' = ' num2str(scan_info(i).Control_PV)], 'fontsize', 20); end;
        SHOT = text(0.23, 1.5, ['Shot #' num2str(j)], 'fontsize', 25);
        
   h_text2 = axes('Position', [0.07, 0.65, 0.3, 0.05], 'Visible', 'off');
        if is_qs_scan; B5D36_text = text(0., 1.7, ['B5D36 = ' num2str(B5D36, '%.2f') ' GeV'], 'fontsize', 20); end;
        nbslice = text(0.04, 0., ['Slices Nb = ' num2str(nb_slice)], 'fontsize',25,'color', 'b');
        if ~is_qs_scan && ~is_qsbend_scan;
            QS_text = text(1.0, 2.5, ['QS = ' num2str(QS, '%.0f') ' GeV'], 'fontsize', 20);
            B5D36_text = text(0., 1.5, ['B5D36 = ' num2str(B5D36, '%.2f') ' GeV'], 'fontsize', 20);
        end;

   h_tex3 = axes('Position', [0.65, 0.1, 0.2, 0.05], 'Visible', 'off');     
        fit = text(0.5,7.,['Expected Fit :'], 'fontsize',18);
        equation = text(0.,5.7,['\sigma_x^2(y) = \sigma_0^2.R_{11}(0)^2 + [4.T_{126}(0)^2.\Theta^2 /\eta_y^2].(y-y_0)^2' ], 'fontsize',15);
        div = text(0.15,4.,['\Theta = ' num2str(1000000*theta2) ' \murad' ], 'fontsize', 25, 'color', 'b');
        ssigma0 = text(0.35,2.7,['\sigma_0 = ' num2str(1000*sigma02) ' \mum' ], 'fontsize', 20, 'color', 'b');
        yy0 = text(0.35,1.5,['y_0 = ' num2str(y02) ' mm' ], 'fontsize', 20, 'color', 'b');
        emittance = text(0.15,0.,['\epsilon < ' num2str(emit2) ' microns' ], 'fontsize', 25, 'color', 'b');
                
        
%%  DISPLAY SETTINGS FOR IMAGES AND PLOTS

LeftAxes1 = axes('position', [0.37, 0.565, 0.22, 0.39]) ;           % CEGAIN display
        image(CEGAIN.yy,1:1392,log10(CEGAIN.img2'),'CDataMapping','scaled');
        colormap(cmap.wbgyr);
        fig_celoss = get(gca,'Children');
        axis xy;
        caxis(CEGAIN_caxis);
        xlabel('x (mm)'); ylabel('y (mm)'); ylim([49 1392]);
        axesPosition = get(gca, 'Position');
        RightAxes1 = axes('Position', axesPosition, 'Color', 'none', 'YAxisLocation', 'right', 'XTick', [], ...
            'Box', 'off','YLim', [49,1392], 'YTick', 50:100:1392, ...
            'YTickLabel', E_EGAIN(50:100:1392));
        ylabel('E (GeV)');
        title('CEGAIN (log scale)');    

LeftAxes2 = axes('position', [0.37, 0.065, 0.22, 0.39]);         % CEGAIN.slice display : desired zoom of the butterfly's parabola
        fig_celoss_slice = image(CEGAIN.yyzoom,xzoom_min:xzoom_max,CEGAIN.slice','CDataMapping','scaled');
        colormap(cmap.wbgyr);
        axis xy;
%         caxis(CEGAIN_caxis);
        xlabel('x (mm)'); ylabel('y (mm)'); ylim([xzoom_min xzoom_max]);
        axesPosition = get(gca, 'Position');
        RightAxes2 = axes('Position', axesPosition, 'Color', 'none', 'YAxisLocation', 'right', 'XTick', [], ...
            'Box', 'off','YLim', [xzoom_min,xzoom_max], 'YTick', xzoom_min:50:xzoom_max, ...
            'YTickLabel', E_EGAIN(xzoom_min:50:xzoom_max));
        ylabel('E (GeV)');
        title('CEGAIN BUTTERFLY (lin scale)');
                       
axes('position', [0.065, 0.065, 0.22, 0.39])  % CEGAIN.reduc display : line-out with first, middle, and last slice
        hold all;
        celossreduc3 = plot(CEGAIN.yyzoom,CEGAIN.reduc3);
        celossreduc2 = plot(CEGAIN.yyzoom,CEGAIN.reduc2);
        celossreduc1 = plot(CEGAIN.yyzoom,CEGAIN.reduc1);
        hleg2 = legend('First Slice','Middle Slice', 'Last Slice');
        set(hleg2,'Location','NorthWest')
        xlabel('x (mm)'); ylabel('Number of Counts');
        title('LINE-OUT');
       
axes('position', [0.7, 0.565, 0.22, 0.39])  % Parabola display: CEGAIN.slice calculations and parabolic fit
        hold all;
        parabola = plot(yparabole2, A2, 'x');  
        parabola_fit = plot(x3,y3);
        hleg1 = legend('Calculations','Parabolic Fit');
        set(hleg1,'Location','North')
        xlabel('y (mm)'); ylabel('\sigma_x^2 (mm^2)'); 
        title('PARABOLIC FIT');
     
           
    else
                          
%%     RESET EACH IMAGE AT EACH ITERATION 
                    
        if is_scan; set(STEP, 'String', [char(regexprep(scan_info(i).Control_PV_name, '_', '\\_')) ' = ' num2str(scan_info(i).Control_PV)]); end;
        
        set(SHOT, 'String', ['Shot #' num2str(j)]);
        set(nbslice,'String', ['Slices Nb = ' num2str(nb_slice)]);
        set(div,'String',['\Theta = ' num2str(1000000*theta2) ' \murad' ]);
        set(yy0,'String',['y_0 = ' num2str(y02) ' mm']);
        set(ssigma0,'String',['\sigma_0 = ' num2str(1000*sigma02) ' \mum' ]);
        set(emittance,'String',['\epsilon < ' num2str(emit2) ' microns']);
        
        set(LeftAxes1,'YLim', [49 1392],'YTick', 50:100:1392, 'YTickLabel', CEGAIN.xx(50:100:1392));
        set(fig_celoss,'XData',CEGAIN.yy,'CData',log10(CEGAIN.img2'));
        set(RightAxes1,'YLim', [49 1392], 'YTick', 50:100:1392,'YTickLabel', E_EGAIN(50:100:1392));
        set(LeftAxes2,'YLim', [xzoom_min xzoom_max],'YTick', xzoom_min:50:xzoom_max, 'YTickLabel', CEGAIN.xx(xzoom_min:50:xzoom_max));
        set(fig_celoss_slice,'XData',CEGAIN.yyzoom,'YData',[xzoom_min xzoom_max],'CData',log10(CEGAIN.slice'));
        set(RightAxes2,'YLim', [xzoom_min,xzoom_max],'YTick', xzoom_min:50:xzoom_max, 'YTickLabel', E_EGAIN(xzoom_min:50:xzoom_max));
        set(celossreduc3,'XData',CEGAIN.yyzoom,'YData',CEGAIN.reduc3);
        set(celossreduc2,'XData',CEGAIN.yyzoom,'YData',CEGAIN.reduc2);
        set(celossreduc1,'XData',CEGAIN.yyzoom,'YData',CEGAIN.reduc1);
        set(parabola,'XData',yparabole2,'YData',A2);
        set(parabola_fit,'XData',x3,'YData',y3);
        
        
    end
    %%   SAVING IN EACH LOOP
     if do_save           
         filename = ['zoom_' num2str(i, '%.2d') '_' num2str(j, '%.3d') '_' num2str(nb_slice) '.png'];
         print('-f1', [save_path data_set '/frames/' filename], '-dpng', '-r100');
     else
         pause(0.1);
     end
  end
% end

 %% SAVING PROCESSED

clear tmp;
if do_save
    if exist([save_path data_set '.mat'], 'file'); tmp = load([save_path data_set]); end;
    tmp.data.processed = processed;
    save([save_path data_set '.mat'], '-struct', 'tmp');
end

