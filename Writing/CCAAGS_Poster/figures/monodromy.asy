import graph3;

size(200,300,keepAspect=false);

material m = opacity(.65) + gray;
pen bluepen = heavyblue + 2.5;
pen mappen = black+2;

currentprojection=orthographic(10,-12,16);
currentlight=(10,10,5);
triple f(pair t) {return (t.x*cos(t.y),t.x*sin(t.y),t.x^(1/2)*sin(t.y/2)+4);}
triple g(pair t) {return (t.x*cos(t.y),t.x*sin(t.y),0);}
triple l(real t) {return (.4*cos(t),.4*sin(t),0);}
triple lift(real t) {return (.4*cos(t),.4*sin(t),sin(t/2)+4);}
triple d(real t) {return (0,0,3-t);}

surface s=surface(f,(0,0),(1,4pi),4,12,Spline);
surface t=surface(g,(0,0),(1,4pi),4,12,Spline);
path3 u1=graph(l,0-pi/3,2pi/3-pi/3+.05);
path3 u2=graph(l,2pi/3-pi/3,4pi/3-pi/3+.05);
path3 u3=graph(l,4pi/3-pi/3,2pi-pi/3);
path3 lift1=graph(lift,0-pi/3,2pi/3-pi/3+.05);
path3 lift2=graph(lift,2pi/3-pi/3,4pi/3-pi/3+.05);
path3 lift3=graph(lift,4pi/3-pi/3,2pi-pi/3);
path3 proj=graph(d,.7,1.8);

dot(point(u1,0),bluepen + 7);
dot(point(lift1,0),bluepen+7);
dot(point(lift3,100),bluepen+7);

draw(s,meshpen=black,surfacepen=m,render(merge=true));
draw(t,meshpen=black,surfacepen=m,render(merge=true));
draw(u1,bluepen,arrow=Arrow3(size=13bp));
draw(u2,bluepen,arrow=Arrow3(size=13bp));
draw(u3,bluepen);
draw(lift1,bluepen,arrow=Arrow3(size=13bp));
draw(lift2,bluepen,arrow=Arrow3(size=13bp));
draw(lift3,bluepen);
draw(proj,mappen,Arrow3(size=8bp));


currentprojection=
  orthographic(camera=(9.9999995073,-12.000000470287501,15.999998857837502),
               up=(-0.002613757,0.003136509,0.01955666),
               target=(0,-1.776357e-15,0),
               zoom=1);

currentprojection=
  orthographic(camera=(9.9999981363315,-11.9999974096695,15.999995956345497),
               up=(-0.002622996,0.003118567,0.01958647),
               target=(0,-1.776357e-15,-1.776357e-15),
               zoom=1);
