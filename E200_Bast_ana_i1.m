% Script to analyze butterfly from E200 2013 data with QS energy = -8 GeV

% Bastien CLAUSSE

% Create: June 16, 2013
% Last edit: June 27, 2013

%% PATH SETTINGS

addpath('C:/Users/Bastien/Desktop/PRE/Matlab/E200_scripts/facet_daq');
addpath('C:/Users/Bastien/Desktop/PRE/Matlab/E200_scripts/bast_ana');

prefix = 'C:/Users/Bastien/Desktop/PRE/Matlab/Volumes/PWFA_4big';
day = '20130430';
data_set = 'E200_10991';
do_save = 0;
save_path = ['C:/Users/Bastien/Desktop/PRE/Matlab/2013_E200_Data_Analysis/Rubidium_density=2.7e17_halfpyro/Zoom_Butterfly_setQSenergy=-8GeV/' day '/'];


%% SCAN INFO FILE INITIALISATION

cmap = custom_cmap();

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

CELOSS = cam_back.CELOSS;

CELOSS.X_RTCL_CTR = 700;
CELOSS.Y_RTCL_CTR = 500;

CELOSS.xx = 1e-3*CELOSS.RESOLUTION * ( (CELOSS.ROI_X-CELOSS.X_RTCL_CTR+1):(CELOSS.ROI_X+CELOSS.ROI_XNP-CELOSS.X_RTCL_CTR) );
CELOSS.yy = 1e-3*CELOSS.RESOLUTION * ( (CELOSS.ROI_Y-CELOSS.Y_RTCL_CTR+1):(CELOSS.ROI_Y+CELOSS.ROI_YNP-CELOSS.Y_RTCL_CTR) );

B5D36 = getB5D36(E200_state);
QS = getQS(E200_state);

E_ELOSS = E200_cher_get_E_axis('20130423', 'CELOSS', 0, 1:1392, 0, B5D36);

order_fit = 2;

%% FIGURE SETTINGS

fig = figure(1);        % positionnement de la figure
set(fig, 'position', [20, 50, 1340, 630]);
set(fig, 'PaperPosition', [0.25, 2.5, 35, 17]);
set(fig, 'color', 'w');
clf();

%% REQUIRED CALCULATIONS FOR STUDY

% for i=1:n_step
    
i=1;

data = load([path mat_filenames{i}]);
if is_qsbend_scan; B5D36 = 20.35 + scan_info(i).Control_PV; end;

[CELOSS.img, ~, CELOSS.pid] = E200_readImages([prefix scan_info(i).CELOSS]);

CELOSS.img = double(CELOSS.img);

for j=1:size(CELOSS.img,3); CELOSS.img(:,:,j) = CELOSS.img(:,:,j) - cam_back.CELOSS.img(:,:); end;

%% ZOOM SETTINGS

    nb_slice = 25;
    [~,ind]=ismember(min((E_ELOSS-(B5D36 + scan_info(i).Control_PV)).^2),(E_ELOSS-(B5D36 + scan_info(i).Control_PV)).^2);
%     [~,ind]=ismember(min((E_ELOSS-12.35).^2),(E_ELOSS-12.35).^2);
    xzoom_min = ind -90;
    xzoom_max = ind +59;
 
    yzoom_min = 300;
    yzoom_max = 675;
 
 
    
    
%%  
 thetamoy = [];
 for j=1:n_shot
% j = 20;
%     for nb_slice=10:10:30


    fprintf('\n%d ', CELOSS.pid(j));
    
    [CELOSS.ana,CELOSS.ana.img] = Ana_CELOSS_img(E_ELOSS, CELOSS.img(:,:,j));
        
    CELOSS.ana.img(CELOSS.ana.img<1) = 1;
    CELOSS.img2 = CELOSS.img(:,:,j);
    CELOSS.img2(CELOSS.img2<1) = 1;


%%  SLICES AND SIGMA  

    CELOSS.yyzoom = CELOSS.yy(yzoom_min:yzoom_max);
    CELOSS.xxzoom = CELOSS.xx(xzoom_min:xzoom_max);

    CELOSS.zoom = CELOSS.ana.img(yzoom_min:yzoom_max, xzoom_min:xzoom_max);
    
    CELOSS.slice = slice(CELOSS.zoom, nb_slice);
    
    h = floor((xzoom_max - xzoom_min +1)/nb_slice);
    
    CELOSS.reduc = reduc(CELOSS.slice, nb_slice);

    CELOSS.reduc1 = CELOSS.reduc(:,h+1);
    CELOSS.reduc2 = CELOSS.reduc(:,h*0.5*nb_slice + 1);
    CELOSS.reduc3 = CELOSS.reduc(:,1+h*(nb_slice-1));
     
    yparabole=[];
    for k=0:(nb_slice - 1); yparabole = [yparabole, CELOSS.xxzoom(1,1+k*h)];end;
    
    A = sigmacarre(CELOSS.slice', CELOSS.yyzoom, nb_slice);
   
    p = polyfit(yparabole, A,order_fit);
    x2 = yparabole;
    y2 = polyval(p,x2);

      
    eta = 0.001*61*B5D36/(B5D36 + scan_info(i).Control_PV);   % en m
%  eta = 0.001*61*20.35/(20.35 + scan_info(i).Control_PV); 
    theta = sqrt(0.25*p(order_fit - 1)*(eta^2)/(23.05^2));    % en rad
    y0 = -p(order_fit)/(2*p(order_fit - 1));
    sigma0 =  sqrt(p(order_fit + 1)-p(order_fit - 1)*y0^2)/(5.6);
    emit = 40000*sigma0*1000*theta;   
  
    thetamoy = [thetamoy , 1000000*theta];
%%  LOOPS INITIALISATION 
    if i==1 && j==1       
%% AXES ET FIGURE SETTINGS

   h_text = axes('Position', [0.07, 0.80, 0.2, 0.05], 'Visible', 'off');     % positionnenemt des informations liées au butterfly sur la figure
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
        div = text(0.15,4.,['\Theta = ' num2str(1000000*theta) ' \murad' ], 'fontsize', 25, 'color', 'b');
        ssigma0 = text(0.35,2.7,['\sigma_0 = ' num2str(1000*sigma0) ' \mum' ], 'fontsize', 20, 'color', 'b');
        yy0 = text(0.35,1.5,['y_0 = ' num2str(y0) ' mm' ], 'fontsize', 20, 'color', 'b');
        emittance = text(0.15,0.,['\epsilon < ' num2str(emit) ' microns' ], 'fontsize', 25, 'color', 'b');
                
        
%%  DISPLAY SETTINGS FOR IMAGES

LeftAxes1 = axes('position', [0.37, 0.565, 0.22, 0.39]) ;                                    % affichage de CELOSS
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

LeftAxes2 = axes('position', [0.37, 0.065, 0.22, 0.39]);         % affichage du zoom sur butterfly de CELOSS.slice
        fig_celoss_slice = image(CELOSS.yyzoom,xzoom_min:xzoom_max,CELOSS.slice','CDataMapping','scaled');
        colormap(cmap.wbgyr);
        axis xy;
%         caxis(CELOSS_caxis);
        xlabel('x (mm)'); ylabel('y (mm)'); ylim([xzoom_min xzoom_max]);
        axesPosition = get(gca, 'Position');
        RightAxes2 = axes('Position', axesPosition, 'Color', 'none', 'YAxisLocation', 'right', 'XTick', [], ...
            'Box', 'off','YLim', [xzoom_min,xzoom_max], 'YTick', xzoom_min:50:xzoom_max, ...
            'YTickLabel', E_ELOSS(xzoom_min:50:xzoom_max));
        ylabel('E (GeV)');
        title('CELOSS BUTTERFLY (lin scale)');

% LeftAxes2 = axes('position', [0.37, 0.065, 0.22, 0.39]);         % affichage du zoom sur butterfly de CELOSS.slice
%         fig_celoss_slice_zoom = image(CELOSS.yyzoom,Xzoom_min:Xzoom_max,log10(CELOSS.slicezoom'),'CDataMapping','scaled');
%         colormap(cmap.wbgyr);
%         axis xy;
%         caxis(CELOSS_caxis);
%         xlabel('x (mm)'); ylabel('y (mm)'); ylim([Xzoom_min Xzoom_max]);
% %         set(gca, 'YTick', xzoom_min:50:xzoom_max);
% %         set(gca, 'YTickLabel', CELOSS.xx(xzoom_min:50:xzoom_max));
%         axesPosition = get(gca, 'Position');
%         RightAxes2 = axes('Position', axesPosition, 'Color', 'none', 'YAxisLocation', 'right', 'XTick', [], ...
%             'Box', 'off','YLim', [Xzoom_min,Xzoom_max], 'YTick', Xzoom_min:50:Xzoom_max, ...
%             'YTickLabel', E_ELOSS(Xzoom_min:50:Xzoom_max));
%         ylabel('E (GeV)');
%         title('CELOSS ZOOM BUTTERFLY (log scale)');
                       
axes('position', [0.065, 0.065, 0.22, 0.39])  % affichage des gaussiennes décrites par CELOSS.reduc
        hold all;
        celossreduc3 = plot(CELOSS.yyzoom,CELOSS.reduc3);
        celossreduc2 = plot(CELOSS.yyzoom,CELOSS.reduc2);
        celossreduc1 = plot(CELOSS.yyzoom,CELOSS.reduc1);
        hleg2 = legend('First Slice','Middle Slice', 'Last Slice');
        set(hleg2,'Location','NorthEast')
        xlabel('x (mm)'); ylabel('nombre de comptes');
        title('LINE-OUT');
       
axes('position', [0.7, 0.565, 0.22, 0.39])  % affichage du vecteur A
        hold all;
        parabola = plot(yparabole, A, 'x');  
        parabola_fit = plot(x2,y2);
        hleg1 = legend('Calculations','Parabolic Fit');
        set(hleg1,'Location','North')
        xlabel('y (mm)'); ylabel('\sigma_x^2 (mm^2)'); 
        title('PARABOLIC FIT');
     
           
    else
                           
%%     RESET EACH IMAGE AT EACH ITERATION 
        if is_scan; set(STEP, 'String', [char(regexprep(scan_info(i).Control_PV_name, '_', '\\_')) ' = ' num2str(scan_info(i).Control_PV)]); end;
        
        set(SHOT, 'String', ['Shot #' num2str(j)]);
        set(nbslice,'String', ['Slices Nb = ' num2str(nb_slice)]);
        set(div,'String',['\Theta = ' num2str(1000000*theta) ' \murad' ]);
        set(yy0,'String',['y_0 = ' num2str(y0) ' mm']);
        set(ssigma0,'String',['\sigma_0 = ' num2str(1000*sigma0) ' \mum' ]);
        set(emittance,'String',['\epsilon < ' num2str(emit) ' microns']);
        
        set(LeftAxes1,'YLim', [49 1392],'YTick', 50:100:1392, 'YTickLabel', CELOSS.xx(50:100:1392));
        set(fig_celoss,'CData',log10(CELOSS.img2'));
        set(RightAxes1,'YLim', [49 1392], 'YTick', 50:100:1392,'YTickLabel', E_ELOSS(50:100:1392));
        set(LeftAxes2,'YLim', [xzoom_min xzoom_max],'YTick', xzoom_min:50:xzoom_max, 'YTickLabel', CELOSS.xx(xzoom_min:50:xzoom_max));
        set(fig_celoss_slice,'YData',[xzoom_min xzoom_max],'CData',log10(CELOSS.slice'));
        set(RightAxes2,'YLim', [xzoom_min,xzoom_max],'YTick', xzoom_min:50:xzoom_max, 'YTickLabel', E_ELOSS(xzoom_min:50:xzoom_max));
%         set(LeftAxes2,'YLim', [Xzoom_min Xzoom_max],'YTick', Xzoom_min:50:Xzoom_max, 'YTickLabel', CELOSS.xx(Xzoom_min:50:Xzoom_max));
%         set(fig_celoss_slice_zoom,'YData',[Xzoom_min Xzoom_max],'CData',log10(CELOSS.slicezoom'));
%         set(RightAxes2,'YLim', [Xzoom_min,Xzoom_max],'YTick', Xzoom_min:50:Xzoom_max, 'YTickLabel', E_ELOSS(Xzoom_min:50:Xzoom_max));
        set(celossreduc3,'XData',CELOSS.yyzoom,'YData',CELOSS.reduc3);
        set(celossreduc2,'XData',CELOSS.yyzoom,'YData',CELOSS.reduc2);
        set(celossreduc1,'XData',CELOSS.yyzoom,'YData',CELOSS.reduc1);
        set(parabola,'XData',yparabole,'YData',A);
        set(parabola_fit,'XData',x2,'YData',y2);
        
                    
        
    end
    %%   SAVING IN EACH LOOP
     if do_save           
         filename = ['zoom_' num2str(i, '%.2d') '_' num2str(j, '%.3d') '_' num2str(nb_slice) '.png'];
         print('-f1', [save_path data_set '/frames/' filename], '-dpng', '-r100');
     else
         pause(0.1);
     end
 end
 
%  thetamoy(2)=0;
%  thetamoy(4)=0;
%  thetamoy(5)=0;
%  thetamoy(9)=0;
%  thetamoy(12)=0;
%  thetamoy(14)=0;
%  thetamoy(15)=0;
%  thetamoy(16)=0;
%  thetamoy(17)=0;
%  thetamoy(18)=0;
%  thetamoy(19)=0;
%  thetamoy(20)=0;
%  thetamoy(26)=0;
%  thetamoy(30)=0;
%  thetamoy(39)=0;
%  thetamoy(32)=0;
%  thetamoy(33)=0;
%  thetamoy(34)=0;
%  thetamoy(38)=0;
%  thetamoy(39)=0;
%  thetamoy(40)=0;
%  thetamoy(41)=0;
%  thetamoy(43)=0;
%  thetamoy(45)=0;
%  thetamoy(46)=0;
%  T = 0;
%  for i=1:n_shot
% T = T + thetamoy(i);
%  end;
% 
%  T = T/(n_shot-12);
% %  
% end