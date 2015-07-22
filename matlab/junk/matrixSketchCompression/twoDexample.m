function twoDexample(n)
  A = zeros(n,3);
  for i=1:n
    r = rand;
    if r <0.25
      mu = [1,1];
      sigma = [0.25,0;0,0.25];
      A(i,1:2) = mvnrnd(mu,sigma);
      A(i,3) = 1;
    elseif r >= 0.25 && r <0.5
      mu = [1,-1];
      sigma = [0.25,0;0,0.25];
      A(i,1:2) = mvnrnd(mu,sigma);
      A(i,3) = 2;
    elseif r >= 0.5 && r <0.75
      mu = [-1,-1];
      sigma = [0.25,0;0,0.25];
      A(i,1:2) = mvnrnd(mu,sigma);
      A(i,3)=3;
    elseif r >= 0.75
      mu = [-1,1];
      sigma = [0.25,0;0,0.25];
      A(i,1:2) = mvnrnd(mu,sigma);
      A(i,3) = 4;
    end
  end
  AA = A(:,1:2);
  AA = normr(AA);
  IND = A(:,3);
  A = [AA IND];
  dlmwrite('example_2d.txt',A);
  l = n/100;
  A = A(:,1:2);
  B = createSketch(l,A);
  B(~any(B,2),:) = [];
  idx = kmeans(B,4);
  dlmwrite('example_2d_sketch.txt',[B idx]);
end