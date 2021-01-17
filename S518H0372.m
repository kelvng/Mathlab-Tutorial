hold on
%% Requirement 1:Compute the coordinate x and y of every point the heart shape.
a=struct();
t = -pi:0.1:pi;
a.xh= round(5 + 4 * sin(t).^5,1);
t = -pi:0.1:pi;
%disp(a.xh)
a.yh= round(3 * cos(t) - 1.7 * cos(2*t) - cos(3*t) +1,1) ;
%disp(a.yh)
%% Requirement 2:Compute the coordinate x and y of every point in the snowflake shape. Then,
%% assign your results into a.Pxn and a.Pyn variables. The parameters are described in StudentID.txt
IDsnowflake=2;
n1=2;
d1=5;
I1='F-F-F-F-F-F-F-F';
F1='F---F+F+F+F+F+F+F---F';
scale=8;
alpha0_1=0.00;
alpha_1=45.00;
x0_1=0.00;
y0_1=0.00;
s=I1;
for h=1:n1
    s=strrep(s,'F',F1);  
end
angle=alpha0_1;

a.Pyn=y0_1;
a.Pxn=x0_1;

for h=1:length(s)
    switch s(h)
        case '-'
            angle=angle-alpha_1;
            
        case 'F'
            x=a.Pxn(end) - d1 * cos(angle * pi / 180);
            y=a.Pyn(end) + d1 * sin(angle * pi / 180);
            a.Pxn=[ a.Pxn x ];
            a.Pyn=[ a.Pyn y ];
        case '+'
            angle=angle+alpha_1;
            
        case '|'
            angle=angle+180; 
    end
end
a.Pxn = round(a.Pxn,1)
%disp(a.Pxn);
%disp(a.Pxn);
%% Requirement 3:Compute the coordinate x and y of every point in the tree shape. Then, assign
%% your results into a.Px and a.Py variables. The parameters are described in StudentID.txt

% 0. PART - INPUT VARIABLES
%============================

  rule(1).vorher = 'F';
  rule(1).danach = 'FF';
  rule(2).vorher = 'X';
  rule(2).danach = 'F-[[X]+X]+F[+FX]-Xo';

  n_Rules = length(rule);
  alpha = 22.50; 
  length_F = 16;
  length_G = 16;
  axiom = 'X';
  n_Repeats = 3;

% 1. PART - CALCULATE THE STRING

for i = 1:n_Repeats
    axiomINcells = cellstr(axiom');
    for j = 1:n_Rules
        hit = strfind(axiom, rule(j).vorher);
        if (length(hit) >= 1)
            for k = hit
                axiomINcells{k} = rule(j).danach;
            end
        end
    end
    axiom = [];
    for j = 1:length(axiomINcells)
        axiom = [axiom, axiomINcells{j}];
    end
end

% 2. PART - PLOT
%=================
IDtree=3;
n2=3;
d2=10;
I2='X';
F2='F-[[X]+X]+F[+FX]-Xo';
alpha0_2=90.00;
alpha_2=22.50;
x0_2=-10.00;
y0_2=-10.00;
xT = -10;
yT = -10;
aT = 0;
da = alpha/180*pi; 
stkPtr = 1;


for i = 1:length(axiom)
    cmdT = axiom(i);
    switch cmdT
    case 'F'
        newxT = xT + length_F*cos(aT);
        newyT = yT + length_F*sin(aT);
        line([yT newyT], [xT newxT],'color','green', 'linewidth',2); 
        xT = newxT;
        yT = newyT;
    case 'X'
        newxT = xT + length_G*cos(aT);
        newyT = yT + length_G*sin(aT);
        xT = newxT;
        yT = newyT;
    case '+'  
        aT = aT + da;
    case '-' 
        aT = aT - da;
    case '[' 
        stack(stkPtr).xT = xT ;
        stack(stkPtr).yT = yT ;
        stack(stkPtr).aT = aT ;
        stkPtr = stkPtr +1 ;
        
    case ']' 
        stkPtr = stkPtr -1 ;
        xT = stack(stkPtr).xT ;
        yT = stack(stkPtr).yT ;
        aT = stack(stkPtr).aT ;
    case 'o'
        a.Pfx = stack(stkPtr).xT;
        a.Pfy = stack(stkPtr).yT;
        plot(a.Pfy,a.Pfx,'yh','LineWidth',3);
    end
end

s=I2;

for h=1:n2
    
    s=strrep(s, 'F', F2);
    
end

angle=alpha0_2;
a.Py=y0_2;
a.Px=x0_2;

for h=1:length(s)
    switch s(h)
        
        case 'F'
            x=a.Px(end)-d2*cos(angle*pi/180);
            y=a.Py(end)+d2*sin(angle*pi/180);
            a.Px=[a.Px x];
            a.Py=[a.Py y];
            
        case '|'
            angle=angle+180; 
            
        case '+'
            angle=angle+alpha_2;
            
        case '-'
            angle=angle-alpha_2;      
    end
end
%% Requirement 4:Compute the coordinate x and y of every flower. Then, assign your results
%% into a.Pfx and a.Pfy variables.
angle=alpha0_2;
a.Pfy=[];
a.Pfx=[];
Stack=[];
YStack=[];
XStack=[];
for h=1:length(s)
    switch s(h)        
        case 'o'      
            y=a.Py(end); 
            x=a.Px(end); 
            a.Pfy=[ a.Pfy y ]; 
            a.Pfx=[ a.Pfx x ];
        case '['
            y=a.Py(end);
            x=a.Px(end);
            Stack=[Stack angle];  
            YStack=[ YStack y ]; 
            XStack=[ XStack x ];
        case ']' 
            n=length(a.Px);
            nStack=length(XStack);
            
            while(a.Px(n)==XStack(nStack) && a.Py(n)==YStack(nStack))
                n=n-1;
                y=a.Py(n);
                x=a.Px(n);
                a.Py=[a.Py y];
                a.Px=[a.Px x];
                
            end
            angle=Stack(end);
            Stack(end)=[];
            YStack(end)=[]; 
            XStack(end)=[];
            a.Pfx=round(a.Pfx,1);
            a.Pfy=round(a.Pfy,1);
    end
end


%% Requirement 5:Compute the number of the snowflake. The equation is defined as
a.nsnowflake=14;
%% Requirement 6:Create a greeting card:
% a)Adding text Happy New Year 2021 to picture. This text is drawn at xt = 50, yt = 0
text(50,0,'Happy New Year 2021','fontsize',14);
% b)Adding all snowslakes to picture. It notes that every snowsflake is reduced to the corresponding
%   scale value (this value is found in StudentID.txt).The new coordinates of every snowsflake are
%   computed as
dx = randi([-50,50]);
dy = randi([-50,50]);
%c) Adding fractal tree to picture
% 0. PART - INPUT VARIABLES
%============================

  rule(1).vorher = 'F';
  rule(1).danach = 'FF';
  rule(2).vorher = 'X';
  rule(2).danach = 'F-[[X]+X]+F[+FX]-Xo';

  n_Rules = length(rule);
  alpha = 22.50; 
  length_F = 16;
  length_G = 16;
  axiom = 'X';
  n_Repeats = 3;

% 1. PART - CALCULATE THE STRING

for i = 1:n_Repeats
    axiomINcells = cellstr(axiom');
    for j = 1:n_Rules
        hit = strfind(axiom, rule(j).vorher);
        if (length(hit) >= 1)
            for k = hit
                axiomINcells{k} = rule(j).danach;
            end
        end
    end
    axiom = [];
    for j = 1:length(axiomINcells)
        axiom = [axiom, axiomINcells{j}];
    end
end

% 2. PART - PLOT
%=================

xT = -10;
yT = -10;
aT = 0;
da = alpha/180*pi; 
stkPtr = 1;


for i = 1:length(axiom)
    cmdT = axiom(i);
    switch cmdT
    case 'F'
        newxT = xT + length_F*cos(aT);
        newyT = yT + length_F*sin(aT);
        line([yT newyT], [xT newxT],'color','green', 'linewidth',2); 
        xT = newxT;
        yT = newyT;
    case 'X'
        newxT = xT + length_G*cos(aT);
        newyT = yT + length_G*sin(aT);
        xT = newxT;
        yT = newyT;
    case '+'  
        aT = aT + da;
    case '-' 
        aT = aT - da;
    case '[' 
        stack(stkPtr).xT = xT ;
        stack(stkPtr).yT = yT ;
        stack(stkPtr).aT = aT ;
        stkPtr = stkPtr +1 ;
        
    case ']' 
        stkPtr = stkPtr -1 ;
        xT = stack(stkPtr).xT ;
        yT = stack(stkPtr).yT ;
        aT = stack(stkPtr).aT ;
    case 'o'
        a.Pfx = stack(stkPtr).xT;
        a.Pfy = stack(stkPtr).yT;
        plot(a.Pfy,a.Pfx,'yh','LineWidth',3);
    end
    a.Pfx=round(a.Pfx,1);
    a.Pfy=round(a.Pfy,1);
end
% d)Adding the heart shape to picture. It note that the position xh and yh of the heart shape are
%   re-computed as xhnew = 2xh + xt ? 20 and yhnew = 2yh, then fill color in the heart shape. Colors
%   are used to draw each shape, described as
%The heart shape is drawn and filled by ’red’ color and line-width = 1.
%  The snowflake is drawn by ’blue’ color and line-width = 1.
% The tree shape is drawn by ’green’ color and line-width = 2.
% The flower is drawn by z color which based on your student ID. The last digit dn in your

syms t;

k = -pi:0.1:pi;
xhn = 40+(2*(4*(sin(t)^5)+5)+(4*(sin(t)^5)+5)-20);
xh = subs(xhn,t,k);

k = -pi:0.1:pi;
yhn = 2*(3*cos(t)-1.7*cos(2*t)-cos(3*t)+1);
yh = subs(yhn,t,k);

fill(xh,yh, 'r','LineWidth',1);

% a.color
n = mod(7,3);
if n==1
    a.color= 'y'
else
    a.color ='r'
end
% a.flower
dn1=1
if dn1 == 0
    a.flower ='p'
elseif dn1 ==1
    a.flower ='h'
elseif dn1 ==2
    a.flower ='*'
elseif dn1 ==3
    a.flower ='o'
elseif dn1 ==4
    a.flower ='+'
elseif dn1 ==5
    a.flower ='v'
elseif dn1 ==6
    a.flower ='x'
elseif dn1 ==7
    a.flower ='d'
elseif dn1 ==8
    a.flower ='s'
else
    a.flower ='>' 
end
%disp(a.flower);
daspect([1,1,1]);
axis([-70,450,-15,350])

saveas(gcf,'pic518H0372.png')
