function [ A ] = sigmacarre( M , x , nb_slice )

[m,n]=size(M);
h = floor(m/nb_slice);
A=[];

for i=0:(nb_slice - 1)
    
    [a,~,~]=mygaussfit(x,M(1+i*h,:));
    A= [A,a*a];
end

end

