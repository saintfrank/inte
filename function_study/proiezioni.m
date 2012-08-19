
m = load('./patt_data.dat') ;

figure(1) ;
plot( m(:, 1) ,m(:, 12) , '+' ) ;
title ("1a Dimension ");
xlabel ("1a ");
ylabel ("e");
axis( [ min( m(:,1) ) , max(m(:,1)), -0.1, max(m(:,12)) ]  );
print("1", "-dpng")

figure(2);
plot( m(:, 2) ,m(:, 12) , '+' ) ;
title ("2a Dimension ");
xlabel ("2a ");
ylabel ("e");
print("2", "-dpng")

figure(3) ;
plot( m(:, 3) ,m(:, 12) , '+' ) ;
title ("3a Dimension ");
xlabel ("3a ");
ylabel ("e");
axis( [ min( m(:,3) ) , max(m(:,3)), -0.1, max(m(:,12)) ]  );
print("3", "-dpng")

figure(4) ;
plot( m(:, 4) ,m(:, 12) , '+' ) ;
title ("4a Dimension ");
xlabel ("4a ");
ylabel ("e");
axis( [ min( m(:,4) ) , max(m(:,4)), -0.1, max(m(:,12)) ]  );
print("4", "-dpng")

figure(5) ;
plot( m(:, 5) ,m(:, 12) , '+' ) ;
title ("5a Dimension ");
xlabel ("5a ");
ylabel ("e");
axis( [ min( m(:,5) ) , max(m(:,5)), -0.1, max(m(:,12)) ]  );
print("5", "-dpng")

figure(6) ;
plot( m(:, 6) ,m(:, 12) , '+' ) ;
title ("6a Dimension ");
xlabel ("6a ");
ylabel ("e");
axis( [ min( m(:,6) ) , max(m(:,6)), -0.1, max(m(:,12)) ]  );
print("6", "-dpng")

figure(7) ;
plot( m(:, 7) ,m(:, 12) , '+' ) ;
title ("7a Dimension ");
xlabel ("7a ");
ylabel ("e");
axis( [ min( m(:,7) ) , max(m(:,7)), -0.1, max(m(:,12)) ]  );
print("7", "-dpng")


figure(8) ;
plot( m(:, 8) ,m(:, 12) , '+' ) ;
title ("8a Dimension ");
xlabel ("8a ");
ylabel ("e");
axis( [ min( m(:,8) ) , max(m(:,8)), -0.1, max(m(:,12)) ]  );
print("8", "-dpng")

figure(9) ;
plot( m(:, 9) ,m(:, 12) , '+' ) ;
title ("9a Dimension ");
xlabel ("9a ");
ylabel ("e");
axis( [ min( m(:,9) ) , max(m(:,9)), -0.1, max(m(:,12)) ]  );
print("9", "-dpng")

figure(10) ;
plot( m(:, 10) ,m(:, 12) , '+' ) ;
title ("10a Dimension ");
xlabel ("10a ");
ylabel ("e");
axis( [ min( m(:,10) ) , max(m(:,10)), -0.1, max(m(:,12)) ]  );
print("10", "-dpng")

