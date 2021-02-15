% This function is used to draw the paths of the first and last cars when an input (solution
% matrix) is given.
function drawpath(s)
n=1;
for i=3:size(s,1)
    if s(i,1)==0
        break
    end
    n = n + 1;
end

first=zeros(n,2); % array for storing locations of the first car
last=zeros(n,2); % array for storing locations of the last car
for i=2:(n+1)
    for j=2:151
        if (abs(s(i,j))<1e-2) && (abs(s(i,j+1))>1e-2) % to locate the last car, 1e-2 as a threshold due to dispersion in Lax-Wendroff scheme
            last(i-1,1)=s(1,j); % storing the locations of the last car
            last(i-1,2)=s(i,1);
        end
        if (s(i,j)~=0) && (s(i,j+1)==0) % to locate the first car
            first(i-1,1)=s(1,j); % storing the locations of the first car
            first(i-1,2)=s(i,1);
        end
    end
end

figure
hold on
plot(first(:,1),first(:,2),last(:,1),last(:,2));
xlabel('x (ft)');
ylabel('t (s)');
legend('First car','Last car');
xlim([-2000 5000]);
ylim([0 80]);
hold off