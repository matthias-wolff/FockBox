%% The Quantum Mouse - Part III
%% Fully Deterministic Perception-Action Cycle
% _<https://www.b-tu.de/en/fg-kommunikationstechnik/team/staff/prof-dr-ing-habil-matthias-wolff 
% Matthias Wolff>_
% 
% <https://www.b-tu.de/en/ BTU Cottbus-Senftenberg>, <https://www.b-tu.de/en/fg-kommunikationstechnik 
% Chair of Communications Engineering>

fock.addFockBoxPaths();
%% 
% 
%% 1   Perception-Action Cycle
% We consider a mouse agent $\mathcal{A}$ which is interconnected with a $N\times 
% M$ maze world $\mathcal{W}$ it lives in by a _perception-action cycle_ shown 
% in Fig. 1.
% 
% 
% 
% *Fig. 1: Perception-action cycle between mouse and maze*
% 
% 
% 
% The mouse can input _actions_
% 
% $$|a\rangle \in\{|\!+\!X\rangle,|0X\rangle,|\!-\!X\rangle\}\otimes\{|\!+\!Y\rangle,|0Y\rangle,|\!-\!Y\rangle\}$$
% 
% into the world. The world, in turn, will output _observations_
% 
% $$|o\rangle\in\mathcal{H}_X\otimes\mathcal{H}_Y\otimes\mathcal{H}_C\otimes\mathcal{H}_G$$
% 
% where $\mathcal{H}_\cdot$ denotes Hilbert spaces
% 
% $\mathcal{H}_X = \mbox{span}\big(|1_X\rangle, \ldots, |N_X\rangle\big)$   
% - x-coordinate space,
% 
% $\mathcal{H}_Y = \mbox{span}\big(|1_Y\rangle, \ldots, |M_Y\rangle\big)$   
% - y-coordinate space,
% 
% $\mathcal{H}_C = \mbox{span}\big(|0_C\rangle, |1_C\rangle\big)$          
% - cheese space, and
% 
% $\mathcal{H}_G = \mbox{span}\big(|0_G\rangle, |1_G\rangle\big)$          
% - water space.
% 
% Throughout this example we assume a fully deterministic stetting, i.e. 
% there is no sensor or actuator noise and no randmoness in the world.
% 
% *1.1   Mouse and Maze*
% 
% First, we create a maze world (Fig. 2):
%%
world = Maze();
world.mmp.setCamPos(-5,-18,23);
for x=1:3
  for y=1:3
    world.addPlace(x,y);
  end
end
for i=1:3
  world.addWall(i,1,'S');
  world.addWall(i,3,'N');
  world.addWall(1,i,'W');
  world.addWall(3,i,'E');
  world.addWater(3,i);
end
world.addWall(2,2,'W');
world.addWall(2,2,'N');
world.addWall(3,2,'N');
world.addCheese(3,3);
title(world.mmp.ax,"Figure 2 - Maze world");
%% 
% Then we create a mouse agent

mouse = Agent.createMouse();
mouse.mmp.setCamPos(-5,-18,23);
title(mouse.mmp.ax,"Figure 2 - Inner stage of mouse agent");
%% 
% and equip it with some (instinctive) knowledge, here four basic movement 
% actions $|a_E\rangle=|{\!+\!}_X\,0_Y\rangle$, $|a_W\rangle=|{\!-\!}_X\,0_Y\rangle$, 
% $|a_N\rangle=|0_X\,{\!+\!}_Y\rangle$, and $|a_N\rangle=|0_X\,{\!-\!}_Y\rangle$:

mouse.addAction('E',fockobj.bket('+X§0Y'));
mouse.addAction('W',fockobj.bket('-X§0Y'));
mouse.addAction('N',fockobj.bket('0X§+Y'));
mouse.addAction('S',fockobj.bket('0X§-Y'));
%% 
% Finally, we set the mouse into the maze at location $(1,1)$

world.addAgent(mouse,1,1);
%% 
% *1.2   Exploring the World*
% 
% [*TODO:* ...]
%%
% Undock figures
set(world.mmp.ax.Parent,'Visible','on','WindowStyle','normal');
set(mouse.mmp.ax.Parent,'Visible','on','WindowStyle','normal');
set(world.mmp.ax.Parent,'OuterPosition',[0,100,500,500]);
set(mouse.mmp.ax.Parent,'OuterPosition',[500,100,500,500]);

%mouse.rexplore(1000);
mouse.sexplore(1:3,1:3);
%% 
% *1.3  The Result of Exploration*
% 
% [*TODO:* ...]
% 
% 
%% 2   Operators on the Inner Stage
% [*TODO:* ...]
%%
PV   = mouse.getVeridicalityProj(); disp(PV);                                   % Veridicality projector
Fxy  = mouse.getFocusOp('XY');      disp(Fxy);                                  % Focus-on-XY operator
UGxy = mouse.getUnfocusOp('XY');    disp(UGxy);                                 % Generic unfocus-from-XY operator
Uxy  = PV*UGxy; Uxy.name='U(XY)';   disp(Uxy);                                  % Veridical unfocus-from-XY operator
PL11 = mouse.getLocationProj(1,1);  disp(PL11);                                 % Projector on location (1,1)