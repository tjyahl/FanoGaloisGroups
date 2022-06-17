import graph3;

size(200,300,keepAspect=false);

material m = opacity(.65) + gray;
pen bluepen = heavyblue + 2.5;
pen redpen = heavyred + 2.5;
pen purplepen = heavymagenta + 2.5;
pen mappen = black + 2;

currentprojection=orthographic(10,-12,16);
currentlight=(10,10,5);
triple g(pair t) {return (t.x*cos(t.y),t.x*sin(t.y),0);}
triple l(real t) {return (.4*cos(t),.4*sin(t),0);}
triple lift(real t) {return (.4*cos(t),.4*sin(t),.3*sin(t/2)+7);}
triple liftadjusted(real t) {return (.4*cos(t)+.01,.4*sin(t),.3*sin(t/2)+7);}
triple d(real t) {return (0,0,3-t);}
triple liftloop(real t) {return (.4*cos(t),.4*sin(t),5.5);}
triple liftloop2(real t) {return (.4*cos(t),.4*sin(t),3.5);}
triple d(real t) {return (0,0,3-t);}

surface t=surface(g,(0,0),(1,4pi),4,12,Spline);
path3 u1=graph(l,0-pi/3,2pi/3-pi/3+.05);
path3 u2=graph(l,2pi/3-pi/3,4pi/3-pi/3+.05);
path3 u3=graph(l,4pi/3-pi/3,2pi-pi/3);
path3 lift1=graph(lift,0-pi/3,2pi/3-pi/3+.05);
path3 lift2=graph(lift,2pi/3-pi/3,4pi/3-pi/3+.05);
path3 lift3=graph(lift,4pi/3-pi/3,2pi-pi/3);
path3 lift4=graph(liftadjusted,2pi-pi/3,8pi/3-pi/3+.05);
path3 lift5=graph(lift,8pi/3-pi/3,10pi/3-pi/3+.05);
path3 lift6=graph(lift,10pi/3-pi/3,4pi-pi/3);
path3 lift7=graph(liftloop,0-pi/3,2pi/3-pi/3+.05);
path3 lift8=graph(liftloop,2pi/3-pi/3,4pi/3-pi/3+.05);
path3 lift9=graph(liftloop,4pi/3-pi/3,2pi-pi/3);
path3 lift10=graph(liftloop2,0-pi/3,2pi/3-pi/3+.05);
path3 lift11=graph(liftloop2,2pi/3-pi/3,4pi/3-pi/3+.05);
path3 lift12=graph(liftloop2,4pi/3-pi/3,2pi-pi/3);
path3 proj=graph(d,.7,1.8);
path3 projj=graph(d,-1.8,-1.2);

dot(point(u1,0),bluepen + 7);
dot(point(lift1,0),bluepen+7);
dot(point(lift7,0),purplepen+7);
dot(point(lift10,0),purplepen+7);
dot(point(lift3,100),redpen+7);
dot(point(projj,0),black+4);
dot(point(projj,50),black+4);
dot(point(projj,100),black+4);

draw(u1,bluepen,arrow=Arrow3(size=13bp));
draw(u2,bluepen,arrow=Arrow3(size=13bp));
draw(u3,bluepen);
draw(lift1,bluepen,arrow=Arrow3(size=13bp));
draw(lift2,bluepen,arrow=Arrow3(size=13bp));
draw(lift3,bluepen);
draw(lift4,redpen,arrow=Arrow3(size=13bp));
draw(lift5,redpen,arrow=Arrow3(size=13bp));
draw(lift6,redpen);
draw(lift7,purplepen,arrow=Arrow3(size=13bp));
draw(lift8,purplepen,arrow=Arrow3(size=13bp));
draw(lift9,purplepen);
draw(lift10,purplepen,arrow=Arrow3(size=13bp));
draw(lift11,purplepen,arrow=Arrow3(size=13bp));
draw(lift12,purplepen);
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

currentprojection=
  orthographic(camera=(9.9999975233595,-11.999996250248504,15.999997592940998),
               up=(-0.0004038058,0.000463107,0.02691652),
               target=(0,-3.552714e-15,-1.776357e-15),
               zoom=1);
