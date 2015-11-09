function A=auxilarymatrix(W)
    D=Dhalfmatrix(W);
    A=D*W*D;
end
