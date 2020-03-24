function [gx,hx] = gx_hx(fy,fx,fyp,fxp);
%This program computes the matrices gx and hx that define the first-order approximation 
%of the DSGE model. That is, if 
%E_t[f(yp,y,xp,x)=0, then the solution is of the form
%xp = h(x,sigma) + sigma * eta * ep
%y = g(x,sigma).
%The first-order approximations to the functions g and h around the point (x,sigma)=(xbar,0), where xbar=h(xbar,0), are:
%h(x,sigma) = xbar + hx (x-xbar) 
%and
%g(x,sigma) = ybar + gx * (x-xbar),
%where ybar=g(xbar,0). 
%Inputs: fy fyp fx fxp
%Outputs: gx hx
%Calls solab.m (by Paul Klein)
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%Date July 17, 2001
A = [-fxp -fyp];
B = [fx fy];

%nfy,nfx,nfyp,nfxp

%a = [-nfxp -nfyp];
%b = [nfx nfy];
%nk=size(nfx,2);
%[gx,hx]=solab(a,b,nk);

[gx,hx]=solab(A,B,size(fx,2));