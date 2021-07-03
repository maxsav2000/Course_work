%% точкинеподвижность
b=3.528;
c=-0.799;
f = @(n,a) n.*exp(-(n./a)+1)
k=100;
set_a = [1:1/k:20];
q = 1:0.01:20
plot(set_a,f(set_a,3),q,q);
%% sdfdf
b=3.528;
c=-0.799;
f = @(n,a) n.*exp((b./n)-c./(n.^2)-a);
k=100;
g = @(x) (b+(b.^2-4.*x.*c).^(1/2))./(2*x); 
set_a = [0.5:1/k:1.836];
set_b = [1.836:1/k:10];
N_u = g(set_a);
plot(set_a,f(N_u,set_a),'g',set_b,f(g(set_b),set_b),'r',1.836,f(g(1.836),1.836),'bo');
legend('Асс. устойчивая','Асс. не устойчивая','Потеря устойчивости');
xlabel('a') 
ylabel('f(N^*)') 
%% устойчивость
k = 500;
set_a = [1/k:1/k:3];
g = @(x) (b+(b.^2-4.*x.*c).^(1/2))./(2*x); 
N_u = g(set_a);
proiz = @(x) 1+2*c./(x.^2)-b./x;
plot(set_a,proiz(N_u));
%%
b=3.528;
c=-0.799;
f = @(n,a) n.*exp((b./n)-c./(n.^2)-a);
n_0 = 1;
k_1 = 2800;
k_2 = 1500;
a_0 = 0.1;
delta = 0.001;
plot(0,0);
hold on
for i = 1:1:k_1
    x = zeros(k_2);
    t = n_0;
    for j = 1:1:k_2
        q = f(t,a_0);
        t = q;
        x(j) = t;
    end
    y = x(k_2-290:k_2);
    A=repmat(a_0,291,1);
    plot(A,y,'ko','MarkerSize',1);
    a_0=a_0+delta;
end
xlabel('a') 
ylabel('N_t') 
%% 2cucle
b=3.528;
c=-0.799;
w=2.76;
f = @(n,a) n.*exp(-(n./a)+1)
k=50;
q = zeros(1,k);
u = zeros(1,k);
t = 1:1:k;
d = 0:1/20:10;
N_0 = 90;
for i=1:k
    u(i) = N_0;
    q(i)=f(N_0,w);
    N_0 = q(i);
end

disp(u);
plot(t,u,'--gs',...
    'LineWidth',2,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('t') 
ylabel('N_t') 
%% 3cucle
b=3.528;
c=-0.799;
f = @(n,a) n.*exp((b./n)-c./(n.^2)-a)
k=100;
set_a = [0.1:1/k:20];
q = 0:0.01:20
rrrr = 0.2;
plot(set_a,f(f(set_a,rrrr),rrrr),q,q,set_a,f(set_a,rrrr));
%% ляпунов
q = 300;
N_0 = 0.05;
b=3.528;
c=-0.799;
f = @(n,a) n.*exp((b./n)-c./(n.^2)-a);
proz = @(n,x) exp(b./n-c./(n.^2)-a)+(2*c/(n.^2)-b./n) .*exp(b./n-c./(n.^2)-a);
res = zeros(1,500);
set_a  = [1/100:1/100:5];
for i = 1:500
    sum = 0;
    u_i = N_0;
    for j=1:q
        sum = sum + log(abs(proz(u_i, set_a(i))));
        u_i = f(u_i,set_a(i));
    end
    sum = sum/q;

    res(i) = sum;
end
plot(set_a,res,'k',[0,5],[0,0]);
xlabel('a') 
ylabel('h_(N_1)') 


