% Script to analyze butterfly from E200 2013 data whatever CELOSS or CEGAIN cameras with an auto-zoom on the parabola's minimum

% Bastien CLAUSSE

% Create: June 16, 2013
% Last edit: July 25, 2013

%%   PATH SETTINGS

addpath('C:/Users/Bastien/Desktop/PRE/Matlab/E200_scripts/facet_daq');
addpath('C:/Users/Bastien/Desktop/PRE/Matlab/E200_scripts/bast_ana');

prefix = 'C:/Users/Bastien/Desktop/PRE/Matlab/Volumes/PWFA_4big';
day = '20130428';
data_set = 'E200_10809';
do_save = 0;
save_path = ['C:/Users/Bastien/Desktop/PRE/Matlab/2013_E200_Data_Analysis/AutoAdjust/' day '/'];


%%   SCAN INFO FILE INITIALISATION

cmap = custom_cmap();

CEGAIN_caxis = [0.8 3.2];
CELOSS_caxis = [0.8 3.2];

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
CELOSS = cam_back.CELOSS;

CEGAIN.X_RTCL_CTR = 700;
CEGAIN.Y_RTCL_CTR = 500;
CELOSS.X_RTCL_CTR = 700;
CELOSS.Y_RTCL_CTR = 500;

CEGAIN.xx = 1e-3*CEGAIN.RESOLUTION * ( (CEGAIN.ROI_X-CEGAIN.X_RTCL_CTR+1):(CEGAIN.ROI_X+CEGAIN.ROI_XNP-CEGAIN.X_RTCL_CTR) );
CEGAIN.yy = 1e-3*CEGAIN.RESOLUTION * ( (CEGAIN.ROI_Y-CEGAIN.Y_RTCL_CTR+1):(CEGAIN.ROI_Y+CEGAIN.ROI_YNP-CEGAIN.Y_RTCL_CTR) );
CELOSS.xx = 1e-3*CELOSS.RESOLUTION * ( (CELOSS.ROI_X-CELOSS.X_RTCL_CTR+1):(CELOSS.ROI_X+CELOSS.ROI_XNP-CELOSS.X_RTCL_CTR) );
CELOSS.yy = 1e-3*CELOSS.RESOLUTION * ( (CELOSS.ROI_Y-CELOSS.Y_RTCL_CTR+1):(CELOSS.ROI_Y+CELOSS.ROI_YNP-CELOSS.Y_RTCL_CTR) );

clear processed;

processed.scalars.DIVERGENCE = zeros(n_step, n_shot);
processed.scalars.SIGMA0 = zeros(n_step, n_shot);
processed.scalars.y0 = zeros(n_step, n_shot);
processed.scalars.EMITTANCE = zeros(n_step, n_shot);

B5D36 = getB5D36(E200_state);
QS = getQS(E200_state);

E_EGAIN = E200_cher_get_E_axis('20130423', 'CEGAIN', 0, 1:1392, 0, B5D36);
E_ELOSS = E200_cher_get_E_axis('20130423', 'CELOSS', 0, 1:1392, 0, B5D36);

order_fit = 2;

%%   FIGURE SETTINGS

fig = figure(1);      
set(fig, 'position', [20, 50, 1340, 630]);
set(fig, 'PaperPosition', [0.25, 2.5, 35, 17]);
set(fig, 'color', 'w');
clf();

%%   REQUIRED CALCULATIONS FOR STUDY

for i=1:n_step
    
% i=2;

data = load([path mat_filenames{i}]);
if is_qsbend_scan; B5D36 = 20.35 + scan_info(i).Control_PV; end;

[CEGAIN.img, ~, CEGAIN.pid] = E200_readImages([prefix scan_info(i).CEGAIN]);
[CELOSS.img, ~, CELOSS.pid] = E200_readImages([prefix scan_info(i).CELOSS]);

CEGAIN.img = double(CEGAIN.img);
CELOSS.img = double(CELOSS.img);

for j=1:size(CEGAIN.img,3); CEGAIN.img(:,:,j) = CEGAIN.img(:,:,j) - cam_back.CEGAIN.img(:,:); end;
for j=1:size(CELOSS.img,3); CELOSS.img(:,:,j) = CELOSS.img(:,:,j) - cam_back.CELOSS.img(:,:); end;

%%   ZOOM INITIALISATION SETTINGS

    nb_slice = 200;

 if scan_info(i).Control_PV < 1
     
    [~,ind]=ismember(min((E_ELOSS-(B5D36 + scan_info(i).Control_PV)).^2),(E_ELOSS-(B5D36 + scan_info(i).Control_PV)).^2);
    xzoom_min = ind - 100;
    xzoom_max = ind + 99;
 
    yzoom_min = 300;
    yzoom_max = 630;
 
 else
     
    [~,ind]=ismember(min((E_EGAIN-(B5D36 + scan_info(i).Control_PV)).^2),(E_EGAIN-(B5D36 + scan_info(i).Control_PV)).^2);
    xzoom_min = ind - 100;
    xzoom_max = ind + 99;
    
    yzoom_min = 300;
    yzoom_max = 630;
 
 end
 
%%    
 for j=1:n_shot
% j = 16;

%%   ANALYSIS OF CEGAIN'S AND CELOSS'S PICTURES

    fprintf('\n%d  \t %d\n', CEGAIN.pid(j), CELOSS.pid(j));
    
    [CEGAIN.ana, CEGAIN.ana.img] = Ana_CEGAIN_img(E_EGAIN, CEGAIN.img(:,:,j));
    [CELOSS.ana,CELOSS.ana.img] = Ana_CELOSS_img(E_ELOSS, CELOSS.img(:,:,j));
        
    CEGAIN.ana.img(CEGAIN.ana.img<1) = 1;
    CEGAIN.img2 = CEGAIN.img(:,:,j);
    CEGAIN.img2(CEGAIN.img2<1) = 1;
    CELOSS.ana.img(CELOSS.ana.img<1) = 1;
    CELOSS.img2 = CELOSS.img(:,:,j);
    CELOSS.img2(CELOSS.img2<1) = 1;


%%   SLICES AND SIGMA  
    
    nb_slicef = 40; % it represents the number of slices, and so the number of points you will use for the fit
                    % keep in mind it's better to conserv (xzoom_maxf - xzoom_minf + 1)/nb_slicef integer
                    % this configuration uses 2 pixels per slices
                    
%%   CELOSS ANALYSIS
    
    CELOSS.yyzoom = CELOSS.yy(yzoom_min:yzoom_max);
    CELOSS.xxzoom = CELOSS.xx(xzoom_min:xzoom_max);
%     CELOSS.zoom = CELOSS.img2(yzoom_min:yzoom_max, xzoom_min:xzoom_max); % Without filters
    CELOSS.zoom = CELOSS.ana.img(yzoom_min:yzoom_max, xzoom_min:xzoom_max); % Using filters    
    CELOSS.slice = slice(CELOSS.zoom, nb_slice);
    
    h = floor((xzoom_max - xzoom_min +1)/nb_slice);
    yparabole=[];
    for k=0:(nb_slice - 1); yparabole = [yparabole, CELOSS.xxzoom(1,1+k*h)];end;
    A = sigmacarre(CELOSS.slice', CELOSS.yyzoom, nb_slice); % A is a vector with the different values of sigma squared with the initial large zoom.
      
    [~,indminA]=ismember(min(A),A);
    Y0 = yparabole(indminA);
    [~,indfinal]=ismember(min((CELOSS.xx - Y0).^2),(CELOSS.xx - Y0).^2);
    xzoom_minf = indfinal -55;              % It defines the window of the final zoom around the parabola's minimum and it can be modified
    xzoom_maxf = indfinal +64;              % if a larger or smaller window is wanted for the butterfly's observation. 
                                            % COMMENT: it's better if (xzoom_maxf - xzoom_minf + 1)/nb_slicef is integer.                                  
    if xzoom_minf < 1   % if the minimum's index is at an extremity of CELOSS image (which is a matrix w/ 1392 rows)
        xzoom_minf = indfinal;
        xzoom_maxf = indfinal +79;        
    end;
    if xzoom_maxf > 1392
        xzoom_minf = indfinal -79;
        xzoom_maxf = indfinal;        
    end;
  
    CELOSS.xxzoomf = CELOSS.xx(xzoom_minf:xzoom_maxf); % it sets the final zoom that will be displayed
    CELOSS.zoomf = CELOSS.ana.img(yzoom_min:yzoom_max, xzoom_minf:xzoom_maxf);
    CELOSS.slicef = slice(CELOSS.zoomf, nb_slicef);
    hf = floor((xzoom_maxf - xzoom_minf + 1)/nb_slicef);
    CELOSS.reduc = reduc(CELOSS.slicef, nb_slicef);    % it defines the line-outs that will be plotted.
    CELOSS.reduc1 = CELOSS.reduc(:,hf+1);
    CELOSS.reduc2 = CELOSS.reduc(:,hf*0.5*nb_slicef + 1);
    CELOSS.reduc3 = CELOSS.reduc(:,1+hf*(nb_slicef-1));
    
    yparabolef=[];
    for k=0:(nb_slicef - 1); yparabolef = [yparabolef, CELOSS.xxzoomf(1,1+k*hf)];end;
    Af = sigmacarre(CELOSS.slicef', CELOSS.yyzoom, nb_slicef); % Af is a vector with the different values of sigma squared with the final zoom aroud the parabola's minimum.
   
    p = polyfit(yparabolef, Af,order_fit);   % Fit calculation
    x2 = yparabolef;
    y2 = polyval(p,x2);
    
%%   PARAMETERS VALUES THANKS TO THE PARABOLIC FIT
    
    eta = 0.001*61*B5D36/(B5D36 + scan_info(i).Control_PV);   % in m
%  eta = 0.001*61*20.35/(20.35 + scan_info(i).Control_PV);  % For example in the dataset E200_10854 the B5D36 is set at 11.35 GeV instead of 20.35 Gev
    theta = sqrt(0.25*p(order_fit - 1)*(eta^2)/(23.05^2));    % in rad
    y0 = -p(order_fit)/(2*p(order_fit - 1)); % in mm
    sigma0 =  sqrt(p(order_fit + 1)-p(order_fit - 1)*y0^2)/(5.6); % in mm
    emit = 40000*sigma0*1000*theta; % in microns    

%%   CEGAIN ANALYSIS   (same thing w/ CEGAIN)
        
    CEGAIN.yyzoom = CEGAIN.yy(yzoom_min:yzoom_max);
    CEGAIN.xxzoom = CEGAIN.xx(xzoom_min:xzoom_max);
    CEGAIN.zoom = CEGAIN.ana.img(yzoom_min:yzoom_max, xzoom_min:xzoom_max);   
    CEGAIN.slice = slice(CEGAIN.zoom, nb_slice);
       
    h = floor((xzoom_max - xzoom_min +1)/nb_slice);
    yparabole2=[];
    for k=0:(nb_slice - 1); yparabole2 = [yparabole2, CEGAIN.xxzoom(1,1+k*h)];end; 
    A2 = sigmacarre(CEGAIN.slice', CEGAIN.yyzoom, nb_slice);    
      
    [~,indminA2]=ismember(min(A2),A2);
    Y02 = yparabole2(indminA2);
    [~,indfinal2]=ismember(min((CEGAIN.xx - Y02).^2),(CEGAIN.xx - Y02).^2);
    xzoom_minf2 = indfinal2 -55;
    xzoom_maxf2 = indfinal2 +64;
    
    if xzoom_minf2 < 1
        xzoom_minf2 = indfinal2;
        xzoom_maxf2 = indfinal2 +79;        
    end;
    if xzoom_maxf2 > 1392
        xzoom_minf2 = indfinal2 -79;
        xzoom_maxf2 = indfinal2;        
    end;
  
    CEGAIN.xxzoomf = CEGAIN.xx(xzoom_minf2:xzoom_maxf2);
    CEGAIN.zoomf = CEGAIN.ana.img(yzoom_min:yzoom_max, xzoom_minf2:xzoom_maxf2);
    CEGAIN.slicef = slice(CEGAIN.zoomf, nb_slicef);
    hf2 = floor((xzoom_maxf2 - xzoom_minf2 + 1)/nb_slicef);
    CEGAIN.reduc = reduc(CEGAIN.slicef, nb_slicef);
    CEGAIN.reduc1 = CEGAIN.reduc(:,hf2+1);
    CEGAIN.reduc2 = CEGAIN.reduc(:,hf2*0.5*nb_slicef + 1);
    CEGAIN.reduc3 = CEGAIN.reduc(:,1+hf2*(nb_slicef-1));
    
    yparabolef2=[];
    for k=0:(nb_slicef - 1); yparabolef2 = [yparabolef2, CEGAIN.xxzoomf(1,1+k*hf2)];end;
    Af2 = sigmacarre(CEGAIN.slicef', CEGAIN.yyzoom, nb_slicef);
   
    p2 = polyfit(yparabolef2, Af2,order_fit);
    x3 = yparabolef2;
    y3 = polyval(p2,x3);
    
    eta2 = 0.001*61*B5D36/(B5D36 + scan_info(i).Control_PV);   % en m
%  eta2 = 0.001*61*20.35/(20.35 + scan_info(i).Control_PV); 
    theta2 = sqrt(0.25*p2(order_fit - 1)*(eta2^2)/(23.05^2));    % in rad
    y02 = -p2(order_fit)/(2*p2(order_fit - 1)); % in mm
    sigma02 =  sqrt(p2(order_fit + 1)-p2(order_fit - 1)*y02^2)/(5.6); % in mm 
    emit2 = 40000*sigma02*1000*theta2; % in microns
    
%% PROCESSED

     if scan_info(i).Control_PV < 1
         
         processed.scalars.DIVERGENCE(i,j) = 1000000*theta;  % in urad
         processed.scalars.SIGMA0(i,j) = 1000*sigma0;  % in um
         processed.scalars.y0(i,j) = y0;  % in mm
         processed.scalars.EMITTANCE(i,j) = emit;  % in microns
         
     else
         
         processed.scalars.DIVERGENCE(i,j) = 1000000*theta2;  % in urad
         processed.scalars.SIGMA0(i,j) = 1000*sigma02;  % in um
         processed.scalars.y0(i,j) = y02;  % in mm
         processed.scalars.EMITTANCE(i,j) = emit2;  % in microns
         
     end
        
%%   LOOPS INITIALISATION 
    if i==1 && j==1       
%%   SETTINGS DISPLAY OF SHOT INFORMATION AND FITTED PARAMETERS 

   h_text = axes('Position', [0.07, 0.80, 0.2, 0.05], 'Visible', 'off');    
        if is_scan; STEP = text(0., 0., [char(regexprep(scan_info(i).Control_PV_name, '_', '\\_')) ' = ' num2str(scan_info(i).Control_PV)], 'fontsize', 20); end;
        SHOT = text(0.23, 1.5, ['Shot #' num2str(j)], 'fontsize', 25);
        
   h_text2 = axes('Position', [0.07, 0.65, 0.3, 0.05], 'Visible', 'off');
        if is_qs_scan; B5D36_text = text(0., 1.7, ['B5D36 = ' num2str(B5D36, '%.2f') ' GeV'], 'fontsize', 20); end;
        nbslice = text(0.04, 0., ['Slices Nb = ' num2str(nb_slicef)], 'fontsize',25,'color', 'b');
        if ~is_qs_scan && ~is_qsbend_scan;
            QS_text = text(1.0, 2.5, ['QS = ' num2str(QS, '%.0f') ' GeV'], 'fontsize', 20);
            B5D36_text = text(0., 1.5, ['B5D36 = ' num2str(B5D36, '%.2f') ' GeV'], 'fontsize', 20);
        end;

   h_tex3 = axes('Position', [0.65, 0.1, 0.2, 0.05], 'Visible', 'off');     
        fit = text(0.5,7.,['Expected Fit :'], 'fontsize',18);
        equation = text(0.,5.7,['\sigma_x^2(y) = \sigma_0^2.R_{11}(0)^2 + [4.T_{126}(0)^2.\Theta^2 /\eta_y^2].(y-y_0)^2' ], 'fontsize',15);
        div = text(0.15,4.1,['\Theta = ' num2str(1000000*theta) ' \murad' ], 'fontsize', 25, 'color', 'b');
        ssigma0 = text(0.35,2.7,['\sigma_0 = ' num2str(1000*sigma0) ' \mum' ], 'fontsize', 20, 'color', 'b');
        yy0 = text(0.35,1.5,['y_0 = ' num2str(y0) ' mm' ], 'fontsize', 20, 'color', 'b');
        emittance = text(0.15,0.,['\epsilon < ' num2str(emit) ' microns' ], 'fontsize', 25, 'color', 'b');
                
        
%%   DISPLAY SETTINGS FOR IMAGES AND PLOTS

LeftAxes1 = axes('position', [0.37, 0.565, 0.22, 0.39]) ;         % CELOSS display
        image(CELOSS.yy,1:1392,log10(CELOSS.img2'),'CDataMapping','scaled');
        colormap(cmap.wbgyr);
        fig_celoss = get(gca,'Children');
        axis xy;
        caxis(CELOSS_caxis);
        xlabel('x (mm)'); ylabel('y (mm)'); ylim([49 1392]);
        axesPosition = get(gca, 'Position');
        RightAxes1 = axes('Position', axesPosition, 'Color', 'none', 'YAxisLocation', 'right', 'XTick', [], ...
            'Box', 'off','YLim', [49,1392], 'YTick', 50:100:1392, ...
            'YTickLabel', E_ELOSS(50:100:1392));
        ylabel('E (GeV)');
        title('CELOSS (log scale)');
        title1=get(gca,'Title');

LeftAxes2 = axes('position', [0.37, 0.065, 0.22, 0.39]);         % CELOSS.slicef display : final zoom around the parabola's minimum 
        fig_celoss_slice = image(CELOSS.yyzoom,xzoom_minf:xzoom_maxf,log10(CELOSS.slicef'),'CDataMapping','scaled');
        colormap(cmap.wbgyr);
        axis xy;
%         caxis(CELOSS_caxis);
        xlabel('x (mm)'); ylabel('y (mm)'); ylim([xzoom_minf xzoom_maxf]);
        axesPosition = get(gca, 'Position');
        RightAxes2 = axes('Position', axesPosition, 'Color', 'none', 'YAxisLocation', 'right', 'XTick', [], ...
            'Box', 'off','YLim', [xzoom_minf,xzoom_maxf], 'YTick', xzoom_minf:10:xzoom_maxf, ...
            'YTickLabel', E_ELOSS(xzoom_minf:10:xzoom_maxf));
        ylabel('E (GeV)');
        title('CELOSS (lin scale)');
        title2=get(gca,'Title');
                            
axes('position', [0.065, 0.065, 0.22, 0.39])  % CELOSS.reduc display : line-out with first, middle, and last slice
        hold all;
        celossreduc3 = plot(CELOSS.yyzoom,CELOSS.reduc3);
        celossreduc2 = plot(CELOSS.yyzoom,CELOSS.reduc2);
        celossreduc1 = plot(CELOSS.yyzoom,CELOSS.reduc1);
        hleg2 = legend('First Slice','Middle Slice', 'Last Slice');
        set(hleg2,'Location','NorthEast')
        xlabel('x (mm)'); ylabel('Number of Counts');
        title('LINE-OUT');
       
axes('position', [0.7, 0.565, 0.22, 0.39])  % Parabola display: CELOSS.slicef calculations and parabolic fit
        hold all;
        parabola = plot(yparabolef, Af, 'x');  
        parabola_fit = plot(x2,y2);
        hleg1 = legend('Calculations','Parabolic Fit');
        set(hleg1,'Location','North')
        xlabel('y (mm)'); ylabel('\sigma_x^2 (mm^2)'); 
        title('PARABOLIC FIT');
     
           
    else
%%   RESET EACH IMAGE AND PLOTS AT EACH ITERATION        
        
            if scan_info(i).Control_PV < 1 % If you need to look the butterfly on the CELOSS camera
                

        if is_scan; set(STEP, 'String', [char(regexprep(scan_info(i).Control_PV_name, '_', '\\_')) ' = ' num2str(scan_info(i).Control_PV)]); end;
        
        set(SHOT, 'String', ['Shot #' num2str(j)]);
        set(nbslice,'String', ['Slices Nb = ' num2str(nb_slicef)]);
        set(div,'String',['\Theta = ' num2str(1000000*theta) ' \murad' ]);
        set(yy0,'String',['y_0 = ' num2str(y0) ' mm']);
        set(ssigma0,'String',['\sigma_0 = ' num2str(1000*sigma0) ' \mum' ]);
        set(emittance,'String',['\epsilon < ' num2str(emit) ' microns']);
        
        set(LeftAxes1,'YLim', [49 1392],'YTick', 50:100:1392, 'YTickLabel', CELOSS.xx(50:100:1392));
        set(fig_celoss,'CData',log10(CELOSS.img2'));
        set(RightAxes1,'YLim', [49 1392], 'YTick', 50:100:1392,'YTickLabel', E_ELOSS(50:100:1392));
        set(LeftAxes2,'YLim', [xzoom_minf xzoom_maxf],'YTick', xzoom_minf:10:xzoom_maxf, 'YTickLabel', CELOSS.xx(xzoom_minf:10:xzoom_maxf));
        set(fig_celoss_slice,'YData',[xzoom_minf xzoom_maxf],'CData',log10(CELOSS.slicef'));
        set(RightAxes2,'YLim', [xzoom_minf,xzoom_maxf],'YTick', xzoom_minf:10:xzoom_maxf, 'YTickLabel', E_ELOSS(xzoom_minf:10:xzoom_maxf));
        set(celossreduc3,'XData',CELOSS.yyzoom,'YData',CELOSS.reduc3);
        set(celossreduc2,'XData',CELOSS.yyzoom,'YData',CELOSS.reduc2);
        set(celossreduc1,'XData',CELOSS.yyzoom,'YData',CELOSS.reduc1);
        set(parabola,'XData',yparabolef,'YData',Af);
        set(parabola_fit,'XData',x2,'YData',y2);
        
            else   % If you need to look the butterfly on the CEGAIN camera
                    
        if is_scan; set(STEP, 'String', [char(regexprep(scan_info(i).Control_PV_name, '_', '\\_')) ' = ' num2str(scan_info(i).Control_PV)]); end;
        
        set(SHOT, 'String', ['Shot #' num2str(j)]);
        set(nbslice,'String', ['Slices Nb = ' num2str(nb_slicef)]);
        set(div,'String',['\Theta = ' num2str(1000000*theta2) ' \murad' ]);
        set(yy0,'String',['y_0 = ' num2str(y02) ' mm']);
        set(ssigma0,'String',['\sigma_0 = ' num2str(1000*sigma02) ' \mum' ]);
        set(emittance,'String',['\epsilon < ' num2str(emit2) ' microns']);
        
        set(LeftAxes1,'YLim', [49 1392],'YTick', 50:100:1392, 'YTickLabel', CEGAIN.xx(50:100:1392));
        set(fig_celoss,'XData',CEGAIN.yy,'CData',log10(CEGAIN.img2'));
        set(title1,'String','CEGAIN (log scale)');
        set(RightAxes1,'YLim', [49 1392], 'YTick', 50:100:1392,'YTickLabel', E_EGAIN(50:100:1392));
        set(LeftAxes2,'YLim', [xzoom_minf2 xzoom_maxf2],'YTick', xzoom_minf2:10:xzoom_maxf2, 'YTickLabel', CEGAIN.xx(xzoom_minf2:10:xzoom_maxf2));
        set(fig_celoss_slice,'XData',CEGAIN.yyzoom,'YData',[xzoom_minf2 xzoom_maxf2],'CData',log10(CEGAIN.slicef'));
        set(title2,'String','CEGAIN (lin scale)');
        set(RightAxes2,'YLim', [xzoom_minf2,xzoom_maxf2],'YTick', xzoom_minf2:10:xzoom_maxf2, 'YTickLabel', E_EGAIN(xzoom_minf2:10:xzoom_maxf2));
        set(celossreduc3,'XData',CEGAIN.yyzoom,'YData',CEGAIN.reduc3);
        set(celossreduc2,'XData',CEGAIN.yyzoom,'YData',CEGAIN.reduc2);
        set(celossreduc1,'XData',CEGAIN.yyzoom,'YData',CEGAIN.reduc1);
        set(parabola,'XData',yparabolef2,'YData',Af2);
        set(parabola_fit,'XData',x3,'YData',y3);
                
            end
            
%%   SAVING PICTURES IN EACH LOOP
     if do_save           
         filename = ['autozoom_' num2str(i, '%.2d') '_' num2str(j, '%.3d') '_' num2str(nb_slicef) '.png'];
         print('-f1', [save_path data_set '/frames/' filename], '-dpng', '-r100');
     else
         pause(0.1);
     end
    end   
 end
end

%% SAVING PROCESSED

clear tmp;
if do_save
    if exist([save_path data_set '.mat'], 'file'); tmp = load([save_path data_set]); end;
    tmp.data.processed = processed;
    save([save_path data_set '.mat'], '-struct', 'tmp');
end







