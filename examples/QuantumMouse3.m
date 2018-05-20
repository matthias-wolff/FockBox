addpath('..');
fock.addFockBoxPaths();

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

% Create mouse and equip it with some konwledge (instinct)                      % -------------------------------------
mouse = Agent.createMouse();                                                    % Create mouse
mouse.mmp.setCamPos(-5,-18,23);
mouse.addAction('E',fockobj.bket('+X§0Y'));                                     % Add east action
mouse.addAction('W',fockobj.bket('-X§0Y'));                                     % Add west action
mouse.addAction('N',fockobj.bket('0X§+Y'));                                     % Add north action
mouse.addAction('S',fockobj.bket('0X§-Y'));                                     % Add south action
world.addAgent(mouse,1,1);                                                      % Add mouse to maze

%mouse.rexplore(1000);
mouse.sexplore(1:3,1:3);

% EOF