OPENQASM 2.0;
include "qelib1.inc";
qreg q[70];
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
h q[35];
h q[36];
h q[37];
h q[38];
h q[39];
h q[40];
h q[41];
h q[42];
h q[43];
h q[44];
h q[45];
h q[46];
h q[47];
h q[48];
h q[49];
h q[50];
h q[51];
h q[52];
h q[53];
h q[54];
h q[55];
h q[56];
h q[57];
h q[58];
h q[59];
h q[60];
h q[61];
h q[62];
h q[63];
h q[64];
h q[65];
h q[66];
h q[67];
h q[68];
h q[69];
h q[14];
t q[5];
s q[61];
h q[68];
cx q[68], q[18];
tdg q[18];
cx q[43], q[18];
t q[18];
cx q[68], q[18];
tdg q[18];
cx q[43], q[18];
t q[18];
cx q[43], q[68];
tdg q[68];
cx q[43], q[68];
t q[43];
t q[68];
h q[68];
h q[10];
s q[0];
h q[15];
cx q[3], q[9];
h q[23];
cx q[23], q[45];
tdg q[45];
cx q[42], q[45];
t q[45];
cx q[23], q[45];
tdg q[45];
cx q[42], q[45];
t q[45];
cx q[42], q[23];
tdg q[23];
cx q[42], q[23];
t q[42];
t q[23];
h q[23];
cx q[47], q[12];
h q[38];
s q[57];
t q[9];
cx q[21], q[65];
h q[8];
cx q[32], q[30];
s q[8];
cx q[62], q[61];
cx q[50], q[3];
s q[16];
h q[15];
s q[50];
h q[14];
cx q[14], q[63];
tdg q[63];
cx q[50], q[63];
t q[63];
cx q[14], q[63];
tdg q[63];
cx q[50], q[63];
t q[63];
cx q[50], q[14];
tdg q[14];
cx q[50], q[14];
t q[50];
t q[14];
h q[14];
cx q[43], q[25];
s q[48];
t q[58];
h q[33];
cx q[33], q[34];
tdg q[34];
cx q[67], q[34];
t q[34];
cx q[33], q[34];
tdg q[34];
cx q[67], q[34];
t q[34];
cx q[67], q[33];
tdg q[33];
cx q[67], q[33];
t q[67];
t q[33];
h q[33];
h q[19];
t q[41];
cx q[6], q[0];
s q[42];
s q[38];
s q[20];
t q[68];
t q[4];
h q[14];
s q[4];
cx q[31], q[40];
h q[16];
cx q[16], q[38];
tdg q[38];
cx q[39], q[38];
t q[38];
cx q[16], q[38];
tdg q[38];
cx q[39], q[38];
t q[38];
cx q[39], q[16];
tdg q[16];
cx q[39], q[16];
t q[39];
t q[16];
h q[16];
t q[13];
h q[69];
h q[59];
cx q[59], q[38];
tdg q[38];
cx q[17], q[38];
t q[38];
cx q[59], q[38];
tdg q[38];
cx q[17], q[38];
t q[38];
cx q[17], q[59];
tdg q[59];
cx q[17], q[59];
t q[17];
t q[59];
h q[59];
s q[15];
h q[54];
cx q[54], q[11];
tdg q[11];
cx q[66], q[11];
t q[11];
cx q[54], q[11];
tdg q[11];
cx q[66], q[11];
t q[11];
cx q[66], q[54];
tdg q[54];
cx q[66], q[54];
t q[66];
t q[54];
h q[54];
h q[20];
s q[31];
cx q[47], q[41];
s q[33];
cx q[47], q[23];
s q[44];
h q[7];
cx q[7], q[36];
tdg q[36];
cx q[48], q[36];
t q[36];
cx q[7], q[36];
tdg q[36];
cx q[48], q[36];
t q[36];
cx q[48], q[7];
tdg q[7];
cx q[48], q[7];
t q[48];
t q[7];
h q[7];
t q[14];
cx q[45], q[33];
t q[11];
t q[17];
h q[37];
cx q[37], q[20];
tdg q[20];
cx q[56], q[20];
t q[20];
cx q[37], q[20];
tdg q[20];
cx q[56], q[20];
t q[20];
cx q[56], q[37];
tdg q[37];
cx q[56], q[37];
t q[56];
t q[37];
h q[37];
h q[28];
cx q[63], q[33];
cx q[59], q[2];
h q[52];
cx q[52], q[39];
tdg q[39];
cx q[5], q[39];
t q[39];
cx q[52], q[39];
tdg q[39];
cx q[5], q[39];
t q[39];
cx q[5], q[52];
tdg q[52];
cx q[5], q[52];
t q[5];
t q[52];
h q[52];
cx q[67], q[6];
s q[56];
h q[8];
s q[16];
cx q[9], q[3];
h q[0];
h q[16];
cx q[16], q[36];
tdg q[36];
cx q[47], q[36];
t q[36];
cx q[16], q[36];
tdg q[36];
cx q[47], q[36];
t q[36];
cx q[47], q[16];
tdg q[16];
cx q[47], q[16];
t q[47];
t q[16];
h q[16];
t q[12];
h q[53];
h q[59];
cx q[55], q[60];
s q[14];
h q[48];
cx q[48], q[61];
tdg q[61];
cx q[3], q[61];
t q[61];
cx q[48], q[61];
tdg q[61];
cx q[3], q[61];
t q[61];
cx q[3], q[48];
tdg q[48];
cx q[3], q[48];
t q[3];
t q[48];
h q[48];
t q[32];
t q[56];
t q[39];
h q[27];
cx q[11], q[46];
s q[58];
h q[13];
cx q[13], q[43];
tdg q[43];
cx q[2], q[43];
t q[43];
cx q[13], q[43];
tdg q[43];
cx q[2], q[43];
t q[43];
cx q[2], q[13];
tdg q[13];
cx q[2], q[13];
t q[2];
t q[13];
h q[13];
t q[4];
h q[20];
t q[51];
h q[46];
s q[22];
s q[19];
h q[7];
cx q[7], q[67];
tdg q[67];
cx q[40], q[67];
t q[67];
cx q[7], q[67];
tdg q[67];
cx q[40], q[67];
t q[67];
cx q[40], q[7];
tdg q[7];
cx q[40], q[7];
t q[40];
t q[7];
h q[7];
h q[29];
h q[17];
cx q[17], q[47];
tdg q[47];
cx q[29], q[47];
t q[47];
cx q[17], q[47];
tdg q[47];
cx q[29], q[47];
t q[47];
cx q[29], q[17];
tdg q[17];
cx q[29], q[17];
t q[29];
t q[17];
h q[17];
cx q[2], q[54];
h q[15];
h q[23];
cx q[23], q[53];
tdg q[53];
cx q[50], q[53];
t q[53];
cx q[23], q[53];
tdg q[53];
cx q[50], q[53];
t q[53];
cx q[50], q[23];
tdg q[23];
cx q[50], q[23];
t q[50];
t q[23];
h q[23];
s q[31];
t q[18];
s q[31];
s q[67];
cx q[32], q[30];
s q[15];
s q[62];
t q[43];
h q[11];
t q[5];
s q[1];
h q[69];
h q[57];
cx q[57], q[58];
tdg q[58];
cx q[47], q[58];
t q[58];
cx q[57], q[58];
tdg q[58];
cx q[47], q[58];
t q[58];
cx q[47], q[57];
tdg q[57];
cx q[47], q[57];
t q[47];
t q[57];
h q[57];
h q[32];
cx q[32], q[15];
tdg q[15];
cx q[19], q[15];
t q[15];
cx q[32], q[15];
tdg q[15];
cx q[19], q[15];
t q[15];
cx q[19], q[32];
tdg q[32];
cx q[19], q[32];
t q[19];
t q[32];
h q[32];
h q[69];
cx q[69], q[52];
tdg q[52];
cx q[10], q[52];
t q[52];
cx q[69], q[52];
tdg q[52];
cx q[10], q[52];
t q[52];
cx q[10], q[69];
tdg q[69];
cx q[10], q[69];
t q[10];
t q[69];
h q[69];
h q[55];
cx q[55], q[37];
tdg q[37];
cx q[59], q[37];
t q[37];
cx q[55], q[37];
tdg q[37];
cx q[59], q[37];
t q[37];
cx q[59], q[55];
tdg q[55];
cx q[59], q[55];
t q[59];
t q[55];
h q[55];
h q[56];
cx q[56], q[68];
tdg q[68];
cx q[6], q[68];
t q[68];
cx q[56], q[68];
tdg q[68];
cx q[6], q[68];
t q[68];
cx q[6], q[56];
tdg q[56];
cx q[6], q[56];
t q[6];
t q[56];
h q[56];
h q[41];
h q[48];
s q[43];
h q[16];
h q[69];
t q[68];
s q[1];
cx q[8], q[19];
h q[45];
cx q[21], q[25];
h q[55];
h q[12];
cx q[12], q[27];
tdg q[27];
cx q[23], q[27];
t q[27];
cx q[12], q[27];
tdg q[27];
cx q[23], q[27];
t q[27];
cx q[23], q[12];
tdg q[12];
cx q[23], q[12];
t q[23];
t q[12];
h q[12];
t q[4];
h q[55];
cx q[55], q[27];
tdg q[27];
cx q[35], q[27];
t q[27];
cx q[55], q[27];
tdg q[27];
cx q[35], q[27];
t q[27];
cx q[35], q[55];
tdg q[55];
cx q[35], q[55];
t q[35];
t q[55];
h q[55];
s q[68];
s q[36];
h q[9];
h q[49];
cx q[49], q[52];
tdg q[52];
cx q[59], q[52];
t q[52];
cx q[49], q[52];
tdg q[52];
cx q[59], q[52];
t q[52];
cx q[59], q[49];
tdg q[49];
cx q[59], q[49];
t q[59];
t q[49];
h q[49];
s q[43];
h q[9];
h q[24];
s q[20];
t q[0];
s q[56];
s q[28];
h q[55];
h q[53];
cx q[53], q[34];
tdg q[34];
cx q[47], q[34];
t q[34];
cx q[53], q[34];
tdg q[34];
cx q[47], q[34];
t q[34];
cx q[47], q[53];
tdg q[53];
cx q[47], q[53];
t q[47];
t q[53];
h q[53];
t q[7];
h q[24];
cx q[24], q[20];
tdg q[20];
cx q[40], q[20];
t q[20];
cx q[24], q[20];
tdg q[20];
cx q[40], q[20];
t q[20];
cx q[40], q[24];
tdg q[24];
cx q[40], q[24];
t q[40];
t q[24];
h q[24];
cx q[1], q[41];
h q[36];
cx q[36], q[28];
tdg q[28];
cx q[66], q[28];
t q[28];
cx q[36], q[28];
tdg q[28];
cx q[66], q[28];
t q[28];
cx q[66], q[36];
tdg q[36];
cx q[66], q[36];
t q[66];
t q[36];
h q[36];
s q[16];
cx q[54], q[23];
h q[23];
cx q[23], q[63];
tdg q[63];
cx q[20], q[63];
t q[63];
cx q[23], q[63];
tdg q[63];
cx q[20], q[63];
t q[63];
cx q[20], q[23];
tdg q[23];
cx q[20], q[23];
t q[20];
t q[23];
h q[23];
h q[49];
cx q[49], q[39];
tdg q[39];
cx q[11], q[39];
t q[39];
cx q[49], q[39];
tdg q[39];
cx q[11], q[39];
t q[39];
cx q[11], q[49];
tdg q[49];
cx q[11], q[49];
t q[11];
t q[49];
h q[49];
cx q[49], q[42];
s q[44];
t q[36];
cx q[49], q[29];
cx q[60], q[68];
t q[41];
h q[42];
cx q[42], q[36];
tdg q[36];
cx q[23], q[36];
t q[36];
cx q[42], q[36];
tdg q[36];
cx q[23], q[36];
t q[36];
cx q[23], q[42];
tdg q[42];
cx q[23], q[42];
t q[23];
t q[42];
h q[42];
s q[31];
t q[61];
s q[38];
h q[33];
cx q[33], q[22];
tdg q[22];
cx q[13], q[22];
t q[22];
cx q[33], q[22];
tdg q[22];
cx q[13], q[22];
t q[22];
cx q[13], q[33];
tdg q[33];
cx q[13], q[33];
t q[13];
t q[33];
h q[33];
cx q[23], q[18];
h q[4];
cx q[4], q[9];
tdg q[9];
cx q[3], q[9];
t q[9];
cx q[4], q[9];
tdg q[9];
cx q[3], q[9];
t q[9];
cx q[3], q[4];
tdg q[4];
cx q[3], q[4];
t q[3];
t q[4];
h q[4];
t q[53];
h q[61];
h q[11];
cx q[11], q[47];
tdg q[47];
cx q[64], q[47];
t q[47];
cx q[11], q[47];
tdg q[47];
cx q[64], q[47];
t q[47];
cx q[64], q[11];
tdg q[11];
cx q[64], q[11];
t q[64];
t q[11];
h q[11];
cx q[58], q[25];
h q[57];
cx q[57], q[1];
tdg q[1];
cx q[68], q[1];
t q[1];
cx q[57], q[1];
tdg q[1];
cx q[68], q[1];
t q[1];
cx q[68], q[57];
tdg q[57];
cx q[68], q[57];
t q[68];
t q[57];
h q[57];
s q[21];
h q[19];
cx q[19], q[48];
tdg q[48];
cx q[42], q[48];
t q[48];
cx q[19], q[48];
tdg q[48];
cx q[42], q[48];
t q[48];
cx q[42], q[19];
tdg q[19];
cx q[42], q[19];
t q[42];
t q[19];
h q[19];
h q[57];
t q[66];
h q[63];
s q[1];
h q[39];
cx q[39], q[25];
tdg q[25];
cx q[27], q[25];
t q[25];
cx q[39], q[25];
tdg q[25];
cx q[27], q[25];
t q[25];
cx q[27], q[39];
tdg q[39];
cx q[27], q[39];
t q[27];
t q[39];
h q[39];
cx q[36], q[50];
cx q[15], q[32];
s q[62];
h q[38];
cx q[38], q[64];
tdg q[64];
cx q[59], q[64];
t q[64];
cx q[38], q[64];
tdg q[64];
cx q[59], q[64];
t q[64];
cx q[59], q[38];
tdg q[38];
cx q[59], q[38];
t q[59];
t q[38];
h q[38];
h q[8];
s q[31];
h q[43];
cx q[30], q[16];
cx q[41], q[45];
cx q[15], q[24];
t q[28];
t q[34];
h q[20];
h q[67];
h q[20];
cx q[4], q[23];
t q[58];
h q[0];
h q[54];
cx q[25], q[40];
t q[46];
s q[17];
s q[40];
t q[55];
s q[25];
s q[60];
h q[69];
h q[28];
cx q[28], q[40];
tdg q[40];
cx q[47], q[40];
t q[40];
cx q[28], q[40];
tdg q[40];
cx q[47], q[40];
t q[40];
cx q[47], q[28];
tdg q[28];
cx q[47], q[28];
t q[47];
t q[28];
h q[28];
cx q[18], q[7];
t q[50];
cx q[14], q[56];
s q[48];
h q[53];
h q[15];
cx q[15], q[59];
tdg q[59];
cx q[60], q[59];
t q[59];
cx q[15], q[59];
tdg q[59];
cx q[60], q[59];
t q[59];
cx q[60], q[15];
tdg q[15];
cx q[60], q[15];
t q[60];
t q[15];
h q[15];
t q[55];
h q[47];
h q[60];
t q[11];
s q[2];
cx q[35], q[45];
t q[23];
cx q[0], q[1];
t q[1];
cx q[0], q[1];
x q[0];
cx q[3], q[4];
t q[4];
s q[3];
cx q[3], q[4];
s q[8];
t q[10];
cx q[9], q[10];
x q[11];
s q[14];
cx q[13], q[14];
y q[15];
y q[16];
t q[19];
x q[19];
cx q[20], q[19];
s q[21];
t q[22];
cx q[22], q[21];
t q[23];
cx q[24], q[23];
cx q[28], q[29];
t q[29];
cx q[28], q[29];
x q[28];
z q[30];
x q[31];
cx q[32], q[31];
y q[33];
y q[34];
cx q[63], q[40];
h q[51];
s q[41];
h q[60];
cx q[60], q[36];
tdg q[36];
cx q[40], q[36];
t q[36];
cx q[60], q[36];
tdg q[36];
cx q[40], q[36];
t q[36];
cx q[40], q[60];
tdg q[60];
cx q[40], q[60];
t q[40];
t q[60];
h q[60];
h q[47];
cx q[47], q[44];
tdg q[44];
cx q[64], q[44];
t q[44];
cx q[47], q[44];
tdg q[44];
cx q[64], q[44];
t q[44];
cx q[64], q[47];
tdg q[47];
cx q[64], q[47];
t q[64];
t q[47];
h q[47];
t q[41];
cx q[49], q[56];
h q[60];
cx q[42], q[67];
h q[48];
cx q[46], q[47];
h q[46];
cx q[46], q[57];
tdg q[57];
cx q[48], q[57];
t q[57];
cx q[46], q[57];
tdg q[57];
cx q[48], q[57];
t q[57];
cx q[48], q[46];
tdg q[46];
cx q[48], q[46];
t q[48];
t q[46];
h q[46];
h q[52];
h q[40];
cx q[40], q[64];
tdg q[64];
cx q[51], q[64];
t q[64];
cx q[40], q[64];
tdg q[64];
cx q[51], q[64];
t q[64];
cx q[51], q[40];
tdg q[40];
cx q[51], q[40];
t q[51];
t q[40];
h q[40];
h q[48];
t q[39];
h q[35];
cx q[35], q[48];
tdg q[48];
cx q[55], q[48];
t q[48];
cx q[35], q[48];
tdg q[48];
cx q[55], q[48];
t q[48];
cx q[55], q[35];
tdg q[35];
cx q[55], q[35];
t q[55];
t q[35];
h q[35];
h q[40];
s q[69];
t q[44];
t q[65];
s q[37];
h q[56];
h q[64];
t q[57];
s q[52];
s q[50];
h q[59];
s q[46];
cx q[69], q[42];
cx q[50], q[63];
t q[54];
cx q[39], q[48];
s q[54];
h q[35];
cx q[35], q[60];
tdg q[60];
cx q[52], q[60];
t q[60];
cx q[35], q[60];
tdg q[60];
cx q[52], q[60];
t q[60];
cx q[52], q[35];
tdg q[35];
cx q[52], q[35];
t q[52];
t q[35];
h q[35];
