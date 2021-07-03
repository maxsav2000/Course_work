%% биомат 1
mu = 1;
eta = 0.9;
gamma = 3;
odefun = @(t, x) [x(1).*(-mu .* log(x(1))-gamma.*x(2)); x(2).*( -1+x(1)./(eta+x(2)))];
f1 = @(x, y) x.*(-mu .* log(x)-gamma.*y);
f2 = @(x, y) y.*( -1+x./(eta+y));



x0_sz = 20;
y0_sz = 1;
x0 = linspace(0, 1, x0_sz);
y0 = [1];
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [1 20], [x0(i); y0(j)]);
        plot(sol(:, 1), sol(:, 2),'r');
        hold on;
    end
end


x0_sz = 1;
y0_sz = 20;
x0 = [1.2];
y0 = linspace(0, 0.3, y0_sz);
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [0 20], [x0(i); y0(j)]);
        plot(sol(:, 1), sol(:, 2),'r');
        hold on;
    end
end

x0_sz = 1;
y0_sz = 20;
x0 = [0.05];
y0 = linspace(0, 1, y0_sz);
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [0 20], [x0(i); y0(j)]);
        plot(sol(:, 1), sol(:, 2),'r');
        hold on;
    end
end
plot(1,0,'k.','MarkerSize',20);
plot(0,0,'k.','MarkerSize',20);
xlabel('u');
ylabel('v');
axis([0.7 1 -.1 0.2]);
syms x
eqn = exp(-x.*gamma./mu) == x + eta;
V = vpasolve(eqn,x,[0 2]);
U = V+eta;
plot(U,V,'k.','MarkerSize',20)

%% биомат 2
mu = 1;
eta = 0.9;
gamma = 3;
odefun = @(t, x) [x(1).*(-mu .* log(x(1))-gamma.*x(2)); x(2).*( -1+x(1)./(eta+x(2)))];
f1 = @(x, y) x.*(-mu .* log(x)-gamma.*y);
f2 = @(x, y) y.*( -1+x./(eta+y));



x0_sz = 20;
y0_sz = 1;
x0 = linspace(0, 1, x0_sz);
y0 = [1];
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [1 20], [x0(i); y0(j)]);
        quiver(sol(:, 1), sol(:, 2), f1(sol(:, 1), sol(:, 2)), f2(sol(:, 1), sol(:, 2)), 'b', 'MarkerSize',10);
        hold on;
    end
end


x0_sz = 1;
y0_sz = 20;
x0 = [1.2];
y0 = linspace(0, 0.3, y0_sz);
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [0 20], [x0(i); y0(j)]);
        quiver(sol(:, 1), sol(:, 2), f1(sol(:, 1), sol(:, 2)), f2(sol(:, 1), sol(:, 2)), 'b', 'MarkerSize',10);
        hold on;
    end
end

x0_sz = 1;
y0_sz = 20;
x0 = [0.05];
y0 = linspace(0, 1, y0_sz);
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [0 20], [x0(i); y0(j)]);
        quiver(sol(:, 1), sol(:, 2), f1(sol(:, 1), sol(:, 2)), f2(sol(:, 1), sol(:, 2)), 'b', 'MarkerSize',10);
        hold on;
    end
end
plot(1,0,'k.','MarkerSize',20);
plot(0,0,'k.','MarkerSize',20);
xlabel('u');
ylabel('v');
axis([0 1.2 0 1]);

%% биомат3
eta = 0.5;
n=30;
mu = linspace(3,10,n);
gamma = linspace(3,10,n);
hold on;
for i=1:n
    for j=1:n
        syms x
        eqn = exp(-x.*gamma(i)./mu(j)) == x + eta;
        V = vpasolve(eqn,x,[0 2]);
        %plot(V,V+eta,'r.','MarkerSize',20);
        U = V+eta;
        J = [-mu(j).*U-mu(j)-gamma(i)*V, -gamma(i).*U; V./(eta+V), -1+U./((eta+V).^2)];
        [Q,D] = eig(J);
        disp(D);
        if (D(1,1).*D(2,2) > 0 )
            plot(mu(i),gamma(j),'r.','MarkerSize',20);
        end
        if (D(1,1).*D(2,2) < 0 )
            plot(mu(i),gamma(j),'b.','MarkerSize',20);
        end
    end
end
xlabel('\gamma')

ylabel('\mu')
%% Область I
a = 2;
b = 2;

odefun = @(t, x) [x(1)*(x(1)-1)*(a-x(1)) - x(1)*x(2)/(x(1) + 1); x(2)*(-b+x(1)/(1+x(1)))];
f1 = @(x, y) x.*(x-1).*(a-x) - x.*y ./(x+1);
f2 = @(x, y) -b*y + x.*y ./(x+1);


x = linspace(0, 3, 1000);
izo = (x .^ 2 - 1) .* (a - x) .* (x > 1) .* (x < a);
plot(x, izo, 'r', 0, 0, 'o black', a, 0, 'o black');
hold on;

x0_sz = 30;
y0_sz = 1;
x0 = linspace(0, 3, x0_sz);
y0 = [1];
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [0 10], [x0(i); y0(j)]);
        quiver(sol(:, 1), sol(:, 2), f1(sol(:, 1), sol(:, 2)), f2(sol(:, 1), sol(:, 2)), 'b', 'MarkerSize',10);
        hold on;
    end
end

x0_sz = 1;
y0_sz = 20;
x0 = [3];
y0 = linspace(0, 1, y0_sz);
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [0 10], [x0(i); y0(j)]);
        quiver(sol(:, 1), sol(:, 2), f1(sol(:, 1), sol(:, 2)), f2(sol(:, 1), sol(:, 2)), 'b', 'MarkerSize',10);
        hold on;
    end
end

[~, sol] = ode45(odefun, 0:-0.001:-30, [1; 0.00001]);
plot(sol(:, 1), sol(:, 2), 'g', 'LineWidth', 1.5);

axis([0 3 0 1]);
%% Область II
a = 10;
b = 1/(1 + 3/(a+sqrt(a^2+3))) + 0.01;

odefun = @(t, x) [x(1)*(x(1)-1)*(a-x(1)) - x(1)*x(2)/(x(1) + 1); x(2)*(-b+x(1)/(1+x(1)))];
f1 = @(x, y) x.*(x-1).*(a-x) - x.*y ./(x+1);
f2 = @(x, y) -b*y + x.*y ./(x+1);

x_eq = b/(1-b);
y_eq = (x_eq^2 - 1)*(a-x_eq);

h = @(x) (x.^2 - 1).*(a-x);

x = linspace(0, a + 1, 1000);
izo = (x .^ 2 - 1) .* (a - x) .* (x > 1) .* (x < a);
plot(x, izo, 'r', b/(1-b) * ones(1, 100), linspace(0, 180, 100), 'r',  0, 0, 'o black');
hold on;

x0_sz = 10;
y0_sz = 10;
x0 = linspace(x_eq - 2, x_eq + 2, x0_sz);
y0 = linspace(y_eq - 3, y_eq + 3, y0_sz);
for i = 1:x0_sz
    for j = 1:y0_sz
        [~, sol] = ode45(odefun, [0 10], [x0(i); y0(j)]);
         quiver(sol(:, 1), sol(:, 2), f1(sol(:, 1), sol(:, 2)), f2(sol(:, 1), sol(:, 2)), 'b', 'LineWidth',0.1, 'MaxHeadSize', 0.03);
%        plot(sol(:, 1), sol(:, 2));
        hold on;
    end
end
 






: