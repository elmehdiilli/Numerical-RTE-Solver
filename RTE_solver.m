function y=RTE_solver(avg_ang,c,albedo,start,length_y,length_z,step_y,step_z,aperture,K,g,alph,g1,g2,mu,nn,ord)

% K=22;
% g=.93;
% g1=.8838;
% g2=-.9835;
% alph=.9832;
% mu=3.483;
% nn=1.33;
% ord=3;
% c=2.2;
albedo=.83;
b=c*albedo;

% phi=rte_angle_dist_org(alph,g1,g2,mu,nn,1+K/2);
% phi=rte_unif_dist_org(1+K/2);
% start=2;
% step_z=.05;
% length_z=2;
% step_y=.01;
% length_y=2;
% aperture=.1016;

I=(length_y/step_y)+1;

M=floor(aperture/step_y)+1;


r(1)=step_y/2;
s=zeros(1,(M-1)/2+1);
s(1,1)=pi*r(1)^2;
% t_ex=zeros(1,T-1);
for n=2:(M-1)/2+1
r(n)=r(n-1)+step_y;
s(n)=pi*r(n)^2-s(n-1);
end


    
    J=(length_z/step_z)+1;
    q=zeros(I,J,K);

q(((I-1)/2)+1,1,1)=1; %% src at the middle
phi=rte_angle_dist_org_pr_biss(g,alph,g1,g2,mu,nn,1+K/2);

 [tt,w, theta]=weight_biss_pr(phi,K,g,alph,g1,g2,mu,nn,ord);




 [radiance_7]=gauss_rte_biss_p(w,theta,length_y,length_z,step_y,step_z,K,q,b,c);
 
 intensity_7=zeros(I,J);
 
  phi(K+1)=2*pi+0.0015/2;
    for k=1:K
% if ((phi(k)>=0 && phi(k)<pi/2) || (phi(k)>=3*pi/2 && phi(k)<=2*pi+0.0015/2))
if ((phi(k)>=0 && phi(k)<=pi))
intensity_7(:,:)=intensity_7(:,:)+radiance_7(:,:,k)*(phi(k+1)-phi(k));
end
    end


%     fprintf('intensity calculated')


for jj=1:J

    
power_7(jj)=s*intensity_7((I-1)/2+1:(M-1)/2+(I-1)/2+1,jj);
    

end
% power=s*intensity((I+1)/2+1:(M-1)/2+(I+1)/2+1,:);

for ii=1:1:J-start/step_z
    
powertrc_7(ii)=power_7(ii+start/step_z);

end
 tm=toc
z=start:step_z:length_z;


  
semilogy(z,powertrc_7./(pi*10^-6),'r-.')

grid on
beep 


