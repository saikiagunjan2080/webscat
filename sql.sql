

create database webscat;
GRANT ALL ON webscat.* TO webscat@localhost IDENTIFIED BY "webscat";

CREATE table LOGINTABLE(
EMAILADDRESS Varchar(50) PRIMARY KEY,
NAME Varchar(50),
CPASSWORD Varchar(50),
DATEOFBIRTH Varchar(50),
PHONENO Varchar(50));


insert into LOGINTABLE values("saikiagunjan@gmail.com","Gunjan Prasad Saikia","gpsaikia","14/12/1990","9957575372");