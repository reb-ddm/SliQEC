OPENQASM 2.0;
include "qelib1.inc";
qreg q[33];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
h q[2];
cx q[2], q[26];
tdg q[26];
cx q[23], q[26];
t q[26];
cx q[2], q[26];
tdg q[26];
cx q[23], q[26];
t q[26];
cx q[23], q[2];
tdg q[2];
cx q[23], q[2];
t q[23];
t q[2];
h q[2];
h q[15];
t q[21];
h q[29];
h q[14];
cx q[14], q[1];
tdg q[1];
cx q[5], q[1];
t q[1];
cx q[14], q[1];
tdg q[1];
cx q[5], q[1];
t q[1];
cx q[5], q[14];
tdg q[14];
cx q[5], q[14];
t q[5];
t q[14];
h q[14];
t q[1];
h q[24];
cx q[20], q[5];
s q[28];
h q[20];
cx q[20], q[15];
tdg q[15];
cx q[16], q[15];
t q[15];
cx q[20], q[15];
tdg q[15];
cx q[16], q[15];
t q[15];
cx q[16], q[20];
tdg q[20];
cx q[16], q[20];
t q[16];
t q[20];
h q[20];
h q[10];
cx q[10], q[28];
tdg q[28];
cx q[15], q[28];
t q[28];
cx q[10], q[28];
tdg q[28];
cx q[15], q[28];
t q[28];
cx q[15], q[10];
tdg q[10];
cx q[15], q[10];
t q[15];
t q[10];
h q[10];
cx q[3], q[13];
cx q[21], q[24];
h q[24];
cx q[24], q[10];
tdg q[10];
cx q[1], q[10];
t q[10];
cx q[24], q[10];
tdg q[10];
cx q[1], q[10];
t q[10];
cx q[1], q[24];
tdg q[24];
cx q[1], q[24];
t q[1];
t q[24];
h q[24];
t q[10];
h q[18];
cx q[18], q[5];
tdg q[5];
cx q[13], q[5];
t q[5];
cx q[18], q[5];
tdg q[5];
cx q[13], q[5];
t q[5];
cx q[13], q[18];
tdg q[18];
cx q[13], q[18];
t q[13];
t q[18];
h q[18];
t q[9];
h q[25];
h q[9];
cx q[9], q[22];
tdg q[22];
cx q[15], q[22];
t q[22];
cx q[9], q[22];
tdg q[22];
cx q[15], q[22];
t q[22];
cx q[15], q[9];
tdg q[9];
cx q[15], q[9];
t q[15];
t q[9];
h q[9];
t q[5];
t q[0];
s q[26];
h q[0];
h q[13];
s q[1];
cx q[0], q[14];
h q[4];
cx q[23], q[7];
s q[16];
cx q[13], q[10];
cx q[15], q[5];
t q[1];
h q[23];
cx q[23], q[25];
tdg q[25];
cx q[16], q[25];
t q[25];
cx q[23], q[25];
tdg q[25];
cx q[16], q[25];
t q[25];
cx q[16], q[23];
tdg q[23];
cx q[16], q[23];
t q[16];
t q[23];
h q[23];
h q[6];
cx q[5], q[20];
s q[27];
s q[23];
s q[10];
t q[13];
cx q[22], q[12];
s q[13];
h q[4];
s q[22];
s q[3];
s q[20];
h q[3];
cx q[3], q[26];
tdg q[26];
cx q[24], q[26];
t q[26];
cx q[3], q[26];
tdg q[26];
cx q[24], q[26];
t q[26];
cx q[24], q[3];
tdg q[3];
cx q[24], q[3];
t q[24];
t q[3];
h q[3];
h q[8];
cx q[8], q[25];
tdg q[25];
cx q[21], q[25];
t q[25];
cx q[8], q[25];
tdg q[25];
cx q[21], q[25];
t q[25];
cx q[21], q[8];
tdg q[8];
cx q[21], q[8];
t q[21];
t q[8];
h q[8];
s q[24];
h q[4];
cx q[2], q[13];
h q[10];
cx q[10], q[3];
tdg q[3];
cx q[20], q[3];
t q[3];
cx q[10], q[3];
tdg q[3];
cx q[20], q[3];
t q[3];
cx q[20], q[10];
tdg q[10];
cx q[20], q[10];
t q[20];
t q[10];
h q[10];
cx q[11], q[22];
t q[22];
h q[16];
h q[22];
cx q[22], q[11];
tdg q[11];
cx q[24], q[11];
t q[11];
cx q[22], q[11];
tdg q[11];
cx q[24], q[11];
t q[11];
cx q[24], q[22];
tdg q[22];
cx q[24], q[22];
t q[24];
t q[22];
h q[22];
h q[8];
cx q[8], q[21];
tdg q[21];
cx q[29], q[21];
t q[21];
cx q[8], q[21];
tdg q[21];
cx q[29], q[21];
t q[21];
cx q[29], q[8];
tdg q[8];
cx q[29], q[8];
t q[29];
t q[8];
h q[8];
h q[16];
cx q[16], q[23];
tdg q[23];
cx q[8], q[23];
t q[23];
cx q[16], q[23];
tdg q[23];
cx q[8], q[23];
t q[23];
cx q[8], q[16];
tdg q[16];
cx q[8], q[16];
t q[8];
t q[16];
h q[16];
s q[1];
h q[0];
cx q[0], q[18];
tdg q[18];
cx q[14], q[18];
t q[18];
cx q[0], q[18];
tdg q[18];
cx q[14], q[18];
t q[18];
cx q[14], q[0];
tdg q[0];
cx q[14], q[0];
t q[14];
t q[0];
h q[0];
t q[13];
cx q[11], q[22];
s q[22];
h q[1];
s q[15];
s q[23];
cx q[9], q[18];
s q[27];
h q[16];
s q[17];
cx q[6], q[2];
h q[13];
cx q[13], q[10];
tdg q[10];
cx q[29], q[10];
t q[10];
cx q[13], q[10];
tdg q[10];
cx q[29], q[10];
t q[10];
cx q[29], q[13];
tdg q[13];
cx q[29], q[13];
t q[29];
t q[13];
h q[13];
cx q[11], q[14];
t q[28];
t q[19];
h q[29];
h q[15];
h q[11];
cx q[11], q[3];
tdg q[3];
cx q[29], q[3];
t q[3];
cx q[11], q[3];
tdg q[3];
cx q[29], q[3];
t q[3];
cx q[29], q[11];
tdg q[11];
cx q[29], q[11];
t q[29];
t q[11];
h q[11];
h q[16];
t q[19];
s q[23];
cx q[12], q[19];
s q[2];
t q[5];
h q[22];
cx q[22], q[1];
tdg q[1];
cx q[3], q[1];
t q[1];
cx q[22], q[1];
tdg q[1];
cx q[3], q[1];
t q[1];
cx q[3], q[22];
tdg q[22];
cx q[3], q[22];
t q[3];
t q[22];
h q[22];
h q[2];
cx q[2], q[18];
tdg q[18];
cx q[23], q[18];
t q[18];
cx q[2], q[18];
tdg q[18];
cx q[23], q[18];
t q[18];
cx q[23], q[2];
tdg q[2];
cx q[23], q[2];
t q[23];
t q[2];
h q[2];
t q[28];
h q[9];
s q[0];
cx q[24], q[2];
h q[27];
cx q[0], q[1];
s q[1];
cx q[0], q[1];
cx q[2], q[3];
z q[5];
z q[8];
s q[13];
cx q[12], q[13];
h q[27];
s q[19];
h q[28];
cx q[28], q[24];
tdg q[24];
cx q[22], q[24];
t q[24];
cx q[28], q[24];
tdg q[24];
cx q[22], q[24];
t q[24];
cx q[22], q[28];
tdg q[28];
cx q[22], q[28];
t q[22];
t q[28];
h q[28];
t q[23];
t q[20];
cx q[29], q[24];
t q[20];
cx q[16], q[23];
h q[20];
t q[19];
h q[21];
h q[29];
h q[21];
t q[16];
h q[15];
cx q[15], q[28];
tdg q[28];
cx q[17], q[28];
t q[28];
cx q[15], q[28];
tdg q[28];
cx q[17], q[28];
t q[28];
cx q[17], q[15];
tdg q[15];
cx q[17], q[15];
t q[17];
t q[15];
h q[15];
cx q[30], q[5];
cx q[31], q[22];
cx q[32], q[1];
