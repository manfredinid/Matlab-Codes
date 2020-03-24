LENGTH = max(size(v));
if ~exist('HP_LAMBDA'),
    HP_LAMBDA = 6.25;
end;
   HP_mat = [1+HP_LAMBDA, -2*HP_LAMBDA, HP_LAMBDA,              zeros(1,LENGTH-3);
             -2*HP_LAMBDA,1+5*HP_LAMBDA,-4*HP_LAMBDA,HP_LAMBDA, zeros(1,LENGTH-4);
                           zeros(LENGTH-4,LENGTH);
              zeros(1,LENGTH-4),HP_LAMBDA,-4*HP_LAMBDA,1+5*HP_LAMBDA,-2*HP_LAMBDA;     
              zeros(1,LENGTH-3),          HP_LAMBDA,   -2*HP_LAMBDA, 1+HP_LAMBDA  ];
   for iiiii=3:LENGTH-2;
     HP_mat(iiiii,iiiii-2)=HP_LAMBDA;
     HP_mat(iiiii,iiiii-1)=-4*HP_LAMBDA;
     HP_mat(iiiii,iiiii)=1+6*HP_LAMBDA;
     HP_mat(iiiii,iiiii+1)=-4*HP_LAMBDA;
     HP_mat(iiiii,iiiii+2)=HP_LAMBDA;
   end;
vtr=HP_mat\v;
vhp=v-vtr; 

      
      