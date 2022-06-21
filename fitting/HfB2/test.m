natom = 48;
dim = 3;
nd = 2870;
ns = 15;
nd123 = 170;
nd1 = 2;
nd2 = 18;
nd3 = 150;

A = readbin("A.bin");
A = reshape(A, [nd nd]);
b = readbin("b.bin");

Phi2 = readbin("Phi2.bin");
Phi2 = reshape(Phi2, [ns ns]);

Phi3 = readbin("Phi3.bin");
Phi3 = reshape(Phi3, [ns ns]);

gd123 = readbin("gd123.bin");
gd = readbin("gdesc.bin");
gdd = readbin("gddesc.bin");
gdd = reshape(gdd, [dim*natom nd]);

eatom = readbin("eatom.bin");
eatom = reshape(eatom, [natom nd123]);
fatom = readbin("fatom.bin");
fatom = reshape(fatom, [dim*natom nd123]);

fatom2 = readbin("fatom2.bin");
fatom2 = reshape(fatom2, [dim*natom nd2]);
fatom3 = readbin("fatom3.bin");
fatom3 = reshape(fatom3, [dim*natom nd3]);

A0 = readbin("A0.bin");
A0 = reshape(A0, [nd nd]);
b0 = readbin("b0.bin");

Phi20 = readbin("Phi20.bin");
Phi20 = reshape(Phi20, [ns 6]);
Phi30 = readbin("Phi30.bin");
Phi30 = reshape(Phi30, [ns 5]);

gd0 = readbin("gdesc0.bin");
gdd0 = readbin("gddesc0.bin");
gdd0 = reshape(gdd0, [dim*natom nd]);

eatom0 = readbin("eatom0.bin");
eatom0 = reshape(eatom0, [natom nd123]);

fatom20 = readbin("fatom20");
fatom20 = reshape(fatom20, [dim*natom nd2]);
fatom30 = readbin("fatom30");
fatom30 = reshape(fatom30, [dim*natom nd3]);

e = A-A0; max(abs(e(:)))
e = b-b0; max(abs(e(:)))
e=Phi2(:,1:6)-Phi20; max(abs(e(:)))
e=Phi3(:,1:5)-Phi30; max(abs(e(:)))

e = eatom-eatom0; max(abs(e(:)))
e = fatom-gdd0(:,1:nd123); max(abs(e(:)))
e = gd-gd0; max(abs(e(:)))
e = gdd(:,1:181)-gdd0(:,1:181); max(abs(e(:)))

e = fatom-gdd(:,1:nd123); max(abs(e(:)))
e = fatom(:,3:20)-fatom2; max(abs(e(:)))
e = fatom(:,21:end)-fatom3; max(abs(e(:)))

e = fatom2-fatom20; max(abs(e(:)))
e = fatom3-fatom30; max(abs(e(:)))

d2 = gd(3:20);
d3 = gd(21:170);
gd23 = zeros(nd2*nd3,1);
gdd23 = zeros(dim*natom, nd2*nd3);
for i = 1:nd3
    for j = 1:nd2
        m = j + nd2*(i-1);
        gd23(m) = d2(j)*d3(i);
        gdd23(:,m) = d2(j)*fatom3(:,i) + d3(i)*fatom2(:,j);
        [max(abs(gdd23(:,m))) max(abs(gdd(:,170+m)))]
        [gd23(m) gd(170+m)]
        m
        pause
    end
end

d20 = gd0(3:20);
d30 = gd0(21:170);
gd230 = zeros(nd2*nd3,1);
gdd230 = zeros(dim*natom, nd2*nd3);
for i = 1:nd3
    for j = 1:nd2
        m = j + nd2*(i-1);
        gd230(m) = d20(j)*d30(i);
        gdd230(:,m) = d20(j)*fatom30(:,i) + d30(i)*fatom20(:,j);
    end
end
gdd1 = [gdd0(:,1:nd123) gdd230];
e = gdd1-gdd0; max(abs(e(:)))


% for (int m3 = 0; m3<M3; m3++)
%     for (int m2 = 0; m2<M2; m2++)
%     {
%         int m = m2 + M2*m3;
%         d23[m] = d2[m2]*d3[m3];                
%         for (int n=0; n<N; n++)
%             dd23[n + N*m] = d2[m2]*dd3[n + N*m3] + dd2[n + N*m2]*d3[m3];
%     }
% 



