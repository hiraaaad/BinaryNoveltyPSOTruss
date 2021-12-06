function newX=enforce_shape_sym(X,max_bar,opt)
switch max_bar
    case 15
    case 25
        if opt.ExploitSymmetry==1
            X([13 25 8 20])=X([10 22 11 23]);
            X([7 19 14 26])=-X([10 22 11 23]);
            X([16 28 17 29])=-X([10 22 11 23]);
            X([9 15 18])=X(12)*[1 1 1];
        end
    case 47
        if opt.ExploitSymmetry==1
            X(3:4:end)=-X(1:4:end);
            X(4:4:end)=X(2:4:end);
        end
    case 224
        X_const_ind=[ 1 2 3 52*3-[2 1 0]   150 153     (4:16:52)*3-2];
        X_const_val=[ 0 0 10  0 -10 0      0   0       zeros(1,4)];
        X(X_const_ind)=X_const_val;
        X(5:48:150)=X(4:48:150);
        for k=0:3
            N=[2 3 4]+16*k;
            M=[6 5 4]+16*k;
            X(M*3-2)=-X(N*3-2);
            X(M*3-1)=X(N*3-1);
            N=[2:6]+16*k;
            M=[10:-1:6]+16*k;
            X(M*3-2)=-X(N*3-1);
            X(M*3-1)=-X(N*3-2);
            N=[2:9]+16*k;
            M=[10:1:17]+16*k;
            X(M*3-2)=-X(N*3-2);
            X(M*3-1)=-X(N*3-1);
            X(([2:17]+16*k)*3)=X((2+16*k)*3);
        end
%         X(5:48:150)=X(4:48:150);
%         for k=0:3
%            if opt.ExploitSymmetry==1
%                N=[2 3 4]+16*k;
%                M=[6 5 4]+16*k;
%                X(M*3-2)=-X(N*3-2);
%                X(M*3-1)=X(N*3-1);
%                N=[2:6]+16*k;
%                M=[10:-1:6]+16*k;
%                X(M*3-2)=-X(N*3-1);
%                X(M*3-1)=-X(N*3-2);
%                N=[2:9]+16*k;
%                M=[10:1:17]+16*k;
%                X(M*3-2)=-X(N*3-2);
%                X(M*3-1)=-X(N*3-1);
%            end
%            X(([2:17]+16*k)*3)=X((2+16*k)*3);
%         end
end
newX=X;