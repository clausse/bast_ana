function [ IMG_SLICE ] = slice( IMG , nb_slice )
% fonction qui définie et calcul les valeurs des tranches selon y 
m=0;
n=0;

[m,n] = size(IMG);

IMG_SLICE=zeros(m,n);
IMG_SLICE_REDUC=zeros(m,n);
h = floor(n/nb_slice);       % définition du pas



for j=0:(nb_slice-1)        % calcul du vecteur moyen pour chaque tranche j
    
   for k=(1+j*h):(h*(j+1))
    IMG_SLICE(:,1+h*j) = IMG_SLICE(:,1+h*j) + IMG(:,k);
   end
   
   IMG_SLICE(:,1+h*j) = IMG_SLICE(:,1+h*j)/h;
end



for j=0:(nb_slice-1)        % génère les tranches sur la matrice
    for i=1:h
       
        IMG_SLICE(:,i+h*j) = IMG_SLICE(:,1+h*j);
   
    end;
end
    
% IMG_SLICE_REDUC(:,1:h) = IMG_SLICE(:,1:h);
% IMG_SLICE_REDUC(:,(1+h*(0.5*nb_slice-1)):(0.5*h*nb_slice)) = IMG_SLICE(:,(1+h*(0.5*nb_slice-1)):(0.5*h*nb_slice));
% IMG_SLICE_REDUC(:,(1+h*(nb_slice-1)):(h*nb_slice)) = IMG_SLICE(:,(1+h*(nb_slice-1)):(h*nb_slice));

end

