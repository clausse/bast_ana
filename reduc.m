function [ B ] = reduc( A , nb_slice )

[m,n]=size(A);
h = floor(n/nb_slice);
B = zeros(m,n);
   
B(:,(h+1):2*h) = A(:,(h+1):2*h);
B(:,(h*0.5*nb_slice + 1):(0.5*h*nb_slice + h)) = A(:,(h*0.5*nb_slice + 1):(0.5*h*nb_slice + h));
B(:,(1+h*(nb_slice-1)):(h*nb_slice )) = A(:,(1+h*(nb_slice-1)):(h*nb_slice));

end

