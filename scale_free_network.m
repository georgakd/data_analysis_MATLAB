clc;
clear all;
rand('state',1);
N=20;
run=1;
active=[];k=[];
nodes=zeros(N,N);
nodes(1,2)=1;
for i=3:N
    for j=1:i-1
        r=rand();
        p=(sum(nodes(j,:))+sum(nodes(:,j)))/sum(sum(nodes));
        if r<=p
           nodes(i,j)=1; 
        if nodes(i,j)==nodes(j,i)
            nodes(i,j)=0;
        end
        end
    end
    for j=i+1:N
        r=rand();
        p=(sum(nodes(j,:))+sum(nodes(:,j)))/sum(sum(nodes));
        if r<=p
           nodes(i,j)=1; 
        if nodes(i,j)==nodes(j,i)
            nodes(i,j)=0;
        end
        end
    end
end

active(run)=sum(sum(nodes));
fprintf('Number of active connections = %d\n\n',active(run));

for i=1:N
   k(i,run)=sum(nodes(i,:))+sum(nodes(:,i));
%   fprintf('Node %d --> %d connections\n',i,k(i,run));
end

x=[];y=[];
for i=1:N
    x(i)=cos((i-1)*((2*pi)/N));
    y(i)=sin((i-1)*((2*pi)/N));
end

clf;
theta=0:0.1:2*pi;
figure (1);
subplot(2,2,1);

hold on;
set(gca,'DataAspectRatio',[1 1 1]);
plot(x,y,'.r',cos(theta),sin(theta),':k');
%fprintf('\nConnections:\n');
for i=1:N
   for j=1:N
      if nodes(i,j)==1
        % fprintf('%d <--> %d\n',i,j);
         xtemp=[x(i) x(j)];
         ytemp=[y(i) y(j)];
         plot(xtemp,ytemp);
      end
   end
end
hold off;



check=0;
while check~=N
    x=[];y=[];
points=zeros(N,N);
for i=1:N
    x(i)=ceil(N * rand());
    y(i)=ceil(N * rand());
    points(x(i),y(i))=1; 
end
check=sum(sum(points));
end



figure (1);
subplot(2,2,2);

hold on;
plot(x,y,'Or');
         set(gca,'DataAspectRatio',[1 1 1]);
         set(gca,'xtick',1:N);
         set(gca,'ytick',1:N);
         grid on;
for i=1:N
    
   for j=1:N
      if nodes(i,j)==1
         xtemp=[x(i) x(j)];
         ytemp=[y(i) y(j)];
         plot(xtemp,ytemp);

      end
   end
end
hold off;


figure (1);
subplot(2,2,3);

hist(k);


figure (1);
subplot(2,2,4);

[k1,k2]=hist(k);
plot(k2,k1,'o');

%disp(nodes);
