OPENQASM 2.0;
include "qelib1.inc";
qreg q[39];
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
h q[30];
h q[31];
h q[32];
h q[33];
h q[34];
t q[13];
s q[14];
s q[0];
h q[23];
h q[22];
h q[16];
cx q[16], q[20];
tdg q[20];
cx q[7], q[20];
t q[20];
cx q[16], q[20];
tdg q[20];
cx q[7], q[20];
t q[20];
cx q[7], q[16];
tdg q[16];
cx q[7], q[16];
t q[7];
t q[16];
h q[16];
h q[21];
cx q[21], q[34];
tdg q[34];
cx q[7], q[34];
t q[34];
cx q[21], q[34];
tdg q[34];
cx q[7], q[34];
t q[34];
cx q[7], q[21];
tdg q[21];
cx q[7], q[21];
t q[7];
t q[21];
h q[21];
h q[3];
t q[19];
t q[16];
h q[26];
t q[12];
t q[24];
cx q[23], q[19];
h q[11];
cx q[11], q[16];
tdg q[16];
cx q[7], q[16];
t q[16];
cx q[11], q[16];
tdg q[16];
cx q[7], q[16];
t q[16];
cx q[7], q[11];
tdg q[11];
cx q[7], q[11];
t q[7];
t q[11];
h q[11];
h q[12];
cx q[8], q[31];
s q[19];
h q[5];
s q[22];
t q[19];
h q[21];
cx q[21], q[2];
tdg q[2];
cx q[25], q[2];
t q[2];
cx q[21], q[2];
tdg q[2];
cx q[25], q[2];
t q[2];
cx q[25], q[21];
tdg q[21];
cx q[25], q[21];
t q[25];
t q[21];
h q[21];
s q[20];
t q[21];
h q[13];
cx q[13], q[29];
tdg q[29];
cx q[28], q[29];
t q[29];
cx q[13], q[29];
tdg q[29];
cx q[28], q[29];
t q[29];
cx q[28], q[13];
tdg q[13];
cx q[28], q[13];
t q[28];
t q[13];
h q[13];
h q[10];
cx q[10], q[27];
tdg q[27];
cx q[28], q[27];
t q[27];
cx q[10], q[27];
tdg q[27];
cx q[28], q[27];
t q[27];
cx q[28], q[10];
tdg q[10];
cx q[28], q[10];
t q[28];
t q[10];
h q[10];
t q[10];
h q[30];
cx q[34], q[20];
cx q[14], q[15];
h q[26];
h q[33];
cx q[33], q[28];
tdg q[28];
cx q[32], q[28];
t q[28];
cx q[33], q[28];
tdg q[28];
cx q[32], q[28];
t q[28];
cx q[32], q[33];
tdg q[33];
cx q[32], q[33];
t q[32];
t q[33];
h q[33];
h q[15];
cx q[15], q[31];
s q[10];
h q[6];
cx q[6], q[7];
tdg q[7];
cx q[14], q[7];
t q[7];
cx q[6], q[7];
tdg q[7];
cx q[14], q[7];
t q[7];
cx q[14], q[6];
tdg q[6];
cx q[14], q[6];
t q[14];
t q[6];
h q[6];
s q[19];
h q[30];
cx q[30], q[19];
tdg q[19];
cx q[1], q[19];
t q[19];
cx q[30], q[19];
tdg q[19];
cx q[1], q[19];
t q[19];
cx q[1], q[30];
tdg q[30];
cx q[1], q[30];
t q[1];
t q[30];
h q[30];
h q[1];
h q[24];
cx q[24], q[16];
tdg q[16];
cx q[25], q[16];
t q[16];
cx q[24], q[16];
tdg q[16];
cx q[25], q[16];
t q[16];
cx q[25], q[24];
tdg q[24];
cx q[25], q[24];
t q[25];
t q[24];
h q[24];
cx q[34], q[30];
t q[28];
t q[9];
t q[17];
h q[22];
h q[28];
h q[13];
cx q[13], q[27];
tdg q[27];
cx q[3], q[27];
t q[27];
cx q[13], q[27];
tdg q[27];
cx q[3], q[27];
t q[27];
cx q[3], q[13];
tdg q[13];
cx q[3], q[13];
t q[3];
t q[13];
h q[13];
h q[15];
s q[28];
h q[22];
cx q[32], q[19];
t q[32];
h q[22];
cx q[22], q[9];
tdg q[9];
cx q[6], q[9];
t q[9];
cx q[22], q[9];
tdg q[9];
cx q[6], q[9];
t q[9];
cx q[6], q[22];
tdg q[22];
cx q[6], q[22];
t q[6];
t q[22];
h q[22];
s q[0];
h q[24];
s q[14];
h q[24];
t q[9];
s q[6];
s q[12];
h q[22];
t q[27];
h q[1];
cx q[1], q[2];
tdg q[2];
cx q[14], q[2];
t q[2];
cx q[1], q[2];
tdg q[2];
cx q[14], q[2];
t q[2];
cx q[14], q[1];
tdg q[1];
cx q[14], q[1];
t q[14];
t q[1];
h q[1];
h q[11];
h q[21];
h q[14];
cx q[14], q[6];
tdg q[6];
cx q[28], q[6];
t q[6];
cx q[14], q[6];
tdg q[6];
cx q[28], q[6];
t q[6];
cx q[28], q[14];
tdg q[14];
cx q[28], q[14];
t q[28];
t q[14];
h q[14];
t q[26];
cx q[34], q[3];
t q[28];
cx q[6], q[30];
cx q[10], q[24];
h q[21];
cx q[21], q[7];
tdg q[7];
cx q[32], q[7];
t q[7];
cx q[21], q[7];
tdg q[7];
cx q[32], q[7];
t q[7];
cx q[32], q[21];
tdg q[21];
cx q[32], q[21];
t q[32];
t q[21];
h q[21];
h q[4];
cx q[4], q[25];
tdg q[25];
cx q[14], q[25];
t q[25];
cx q[4], q[25];
tdg q[25];
cx q[14], q[25];
t q[25];
cx q[14], q[4];
tdg q[4];
cx q[14], q[4];
t q[14];
t q[4];
h q[4];
t q[11];
t q[26];
cx q[1], q[27];
h q[24];
cx q[20], q[3];
h q[3];
cx q[3], q[8];
tdg q[8];
cx q[1], q[8];
t q[8];
cx q[3], q[8];
tdg q[8];
cx q[1], q[8];
t q[8];
cx q[1], q[3];
tdg q[3];
cx q[1], q[3];
t q[1];
t q[3];
h q[3];
h q[8];
cx q[8], q[5];
tdg q[5];
cx q[4], q[5];
t q[5];
cx q[8], q[5];
tdg q[5];
cx q[4], q[5];
t q[5];
cx q[4], q[8];
tdg q[8];
cx q[4], q[8];
t q[4];
t q[8];
h q[8];
cx q[34], q[33];
s q[0];
cx q[10], q[23];
h q[4];
s q[3];
t q[20];
cx q[30], q[26];
h q[17];
cx q[17], q[10];
tdg q[10];
cx q[19], q[10];
t q[10];
cx q[17], q[10];
tdg q[10];
cx q[19], q[10];
t q[10];
cx q[19], q[17];
tdg q[17];
cx q[19], q[17];
t q[19];
t q[17];
h q[17];
cx q[21], q[12];
h q[4];
cx q[4], q[24];
tdg q[24];
cx q[32], q[24];
t q[24];
cx q[4], q[24];
tdg q[24];
cx q[32], q[24];
t q[24];
cx q[32], q[4];
tdg q[4];
cx q[32], q[4];
t q[32];
t q[4];
h q[4];
cx q[3], q[0];
h q[34];
s q[15];
h q[6];
h q[21];
h q[26];
h q[6];
s q[5];
cx q[12], q[21];
s q[4];
h q[14];
s q[30];
t q[4];
cx q[25], q[32];
cx q[12], q[18];
t q[0];
x q[5];
t q[7];
cx q[6], q[7];
cx q[8], q[9];
s q[8];
z q[10];
y q[13];
cx q[14], q[15];
t q[15];
t q[14];
cx q[14], q[15];
cx q[18], q[25];
t q[24];
h q[19];
cx q[19], q[18];
tdg q[18];
cx q[32], q[18];
t q[18];
cx q[19], q[18];
tdg q[18];
cx q[32], q[18];
t q[18];
cx q[32], q[19];
tdg q[19];
cx q[32], q[19];
t q[32];
t q[19];
h q[19];
h q[18];
h q[25];
s q[30];
s q[29];
h q[30];
s q[23];
cx q[32], q[19];
cx q[31], q[17];
s q[25];
t q[34];
h q[20];
h q[32];
cx q[32], q[30];
tdg q[30];
cx q[18], q[30];
t q[30];
cx q[32], q[30];
tdg q[30];
cx q[18], q[30];
t q[30];
cx q[18], q[32];
tdg q[32];
cx q[18], q[32];
t q[18];
t q[32];
h q[32];
h q[30];
cx q[30], q[21];
tdg q[21];
cx q[31], q[21];
t q[21];
cx q[30], q[21];
tdg q[21];
cx q[31], q[21];
t q[21];
cx q[31], q[30];
tdg q[30];
cx q[31], q[30];
t q[31];
t q[30];
h q[30];
h q[21];
cx q[21], q[22];
tdg q[22];
cx q[25], q[22];
t q[22];
cx q[21], q[22];
tdg q[22];
cx q[25], q[22];
t q[22];
cx q[25], q[21];
tdg q[21];
cx q[25], q[21];
t q[25];
t q[21];
h q[21];
h q[25];
cx q[35], q[26];
cx q[36], q[33];
cx q[37], q[15];
cx q[38], q[6];
