OPENQASM 2.0;
include "qelib1.inc";
qreg q[80];
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
h q[70];
h q[71];
h q[72];
h q[73];
h q[74];
h q[75];
h q[76];
h q[77];
h q[78];
h q[79];
h q[44];
cx q[44], q[66];
s q[63];
cx q[61], q[38];
cx q[75], q[30];
cx q[35], q[3];
h q[77];
cx q[77], q[27];
tdg q[27];
cx q[20], q[27];
t q[27];
cx q[77], q[27];
tdg q[27];
cx q[20], q[27];
t q[27];
cx q[20], q[77];
tdg q[77];
cx q[20], q[77];
t q[20];
t q[77];
h q[77];
cx q[8], q[2];
s q[45];
h q[10];
cx q[10], q[19];
tdg q[19];
cx q[6], q[19];
t q[19];
cx q[10], q[19];
tdg q[19];
cx q[6], q[19];
t q[19];
cx q[6], q[10];
tdg q[10];
cx q[6], q[10];
t q[6];
t q[10];
h q[10];
cx q[79], q[1];
cx q[6], q[21];
t q[75];
h q[37];
cx q[37], q[44];
tdg q[44];
cx q[41], q[44];
t q[44];
cx q[37], q[44];
tdg q[44];
cx q[41], q[44];
t q[44];
cx q[41], q[37];
tdg q[37];
cx q[41], q[37];
t q[41];
t q[37];
h q[37];
t q[26];
t q[8];
t q[24];
t q[48];
s q[41];
h q[29];
cx q[29], q[15];
tdg q[15];
cx q[61], q[15];
t q[15];
cx q[29], q[15];
tdg q[15];
cx q[61], q[15];
t q[15];
cx q[61], q[29];
tdg q[29];
cx q[61], q[29];
t q[61];
t q[29];
h q[29];
t q[5];
s q[1];
h q[0];
cx q[69], q[8];
cx q[36], q[18];
h q[20];
cx q[20], q[78];
tdg q[78];
cx q[60], q[78];
t q[78];
cx q[20], q[78];
tdg q[78];
cx q[60], q[78];
t q[78];
cx q[60], q[20];
tdg q[20];
cx q[60], q[20];
t q[60];
t q[20];
h q[20];
t q[77];
h q[76];
cx q[76], q[50];
tdg q[50];
cx q[19], q[50];
t q[50];
cx q[76], q[50];
tdg q[50];
cx q[19], q[50];
t q[50];
cx q[19], q[76];
tdg q[76];
cx q[19], q[76];
t q[19];
t q[76];
h q[76];
cx q[14], q[58];
s q[31];
t q[39];
t q[64];
t q[32];
t q[59];
t q[4];
h q[11];
cx q[57], q[3];
t q[48];
t q[72];
h q[71];
cx q[71], q[26];
tdg q[26];
cx q[68], q[26];
t q[26];
cx q[71], q[26];
tdg q[26];
cx q[68], q[26];
t q[26];
cx q[68], q[71];
tdg q[71];
cx q[68], q[71];
t q[68];
t q[71];
h q[71];
h q[18];
h q[41];
cx q[41], q[49];
tdg q[49];
cx q[5], q[49];
t q[49];
cx q[41], q[49];
tdg q[49];
cx q[5], q[49];
t q[49];
cx q[5], q[41];
tdg q[41];
cx q[5], q[41];
t q[5];
t q[41];
h q[41];
s q[22];
h q[16];
cx q[16], q[76];
tdg q[76];
cx q[23], q[76];
t q[76];
cx q[16], q[76];
tdg q[76];
cx q[23], q[76];
t q[76];
cx q[23], q[16];
tdg q[16];
cx q[23], q[16];
t q[23];
t q[16];
h q[16];
s q[1];
h q[74];
h q[24];
cx q[24], q[49];
tdg q[49];
cx q[3], q[49];
t q[49];
cx q[24], q[49];
tdg q[49];
cx q[3], q[49];
t q[49];
cx q[3], q[24];
tdg q[24];
cx q[3], q[24];
t q[3];
t q[24];
h q[24];
t q[41];
h q[16];
s q[74];
h q[2];
h q[4];
s q[49];
h q[6];
cx q[6], q[18];
tdg q[18];
cx q[45], q[18];
t q[18];
cx q[6], q[18];
tdg q[18];
cx q[45], q[18];
t q[18];
cx q[45], q[6];
tdg q[6];
cx q[45], q[6];
t q[45];
t q[6];
h q[6];
t q[30];
h q[71];
s q[5];
cx q[79], q[41];
t q[68];
s q[41];
t q[77];
s q[69];
cx q[25], q[63];
h q[10];
cx q[10], q[33];
tdg q[33];
cx q[55], q[33];
t q[33];
cx q[10], q[33];
tdg q[33];
cx q[55], q[33];
t q[33];
cx q[55], q[10];
tdg q[10];
cx q[55], q[10];
t q[55];
t q[10];
h q[10];
h q[42];
h q[7];
t q[25];
s q[78];
cx q[11], q[17];
cx q[7], q[58];
t q[61];
cx q[1], q[34];
s q[61];
t q[36];
t q[7];
cx q[25], q[53];
cx q[40], q[38];
t q[70];
t q[29];
t q[78];
cx q[0], q[38];
h q[3];
s q[68];
cx q[7], q[15];
s q[26];
t q[71];
h q[34];
cx q[34], q[42];
tdg q[42];
cx q[3], q[42];
t q[42];
cx q[34], q[42];
tdg q[42];
cx q[3], q[42];
t q[42];
cx q[3], q[34];
tdg q[34];
cx q[3], q[34];
t q[3];
t q[34];
h q[34];
h q[63];
h q[50];
cx q[50], q[24];
tdg q[24];
cx q[56], q[24];
t q[24];
cx q[50], q[24];
tdg q[24];
cx q[56], q[24];
t q[24];
cx q[56], q[50];
tdg q[50];
cx q[56], q[50];
t q[56];
t q[50];
h q[50];
h q[20];
cx q[20], q[26];
tdg q[26];
cx q[40], q[26];
t q[26];
cx q[20], q[26];
tdg q[26];
cx q[40], q[26];
t q[26];
cx q[40], q[20];
tdg q[20];
cx q[40], q[20];
t q[40];
t q[20];
h q[20];
h q[11];
cx q[11], q[31];
tdg q[31];
cx q[23], q[31];
t q[31];
cx q[11], q[31];
tdg q[31];
cx q[23], q[31];
t q[31];
cx q[23], q[11];
tdg q[11];
cx q[23], q[11];
t q[23];
t q[11];
h q[11];
h q[70];
h q[68];
cx q[68], q[6];
tdg q[6];
cx q[13], q[6];
t q[6];
cx q[68], q[6];
tdg q[6];
cx q[13], q[6];
t q[6];
cx q[13], q[68];
tdg q[68];
cx q[13], q[68];
t q[13];
t q[68];
h q[68];
h q[44];
h q[78];
h q[54];
cx q[54], q[37];
tdg q[37];
cx q[64], q[37];
t q[37];
cx q[54], q[37];
tdg q[37];
cx q[64], q[37];
t q[37];
cx q[64], q[54];
tdg q[54];
cx q[64], q[54];
t q[64];
t q[54];
h q[54];
h q[58];
h q[35];
cx q[35], q[43];
tdg q[43];
cx q[38], q[43];
t q[43];
cx q[35], q[43];
tdg q[43];
cx q[38], q[43];
t q[43];
cx q[38], q[35];
tdg q[35];
cx q[38], q[35];
t q[38];
t q[35];
h q[35];
h q[12];
cx q[12], q[39];
tdg q[39];
cx q[76], q[39];
t q[39];
cx q[12], q[39];
tdg q[39];
cx q[76], q[39];
t q[39];
cx q[76], q[12];
tdg q[12];
cx q[76], q[12];
t q[76];
t q[12];
h q[12];
s q[76];
h q[4];
cx q[4], q[48];
tdg q[48];
cx q[49], q[48];
t q[48];
cx q[4], q[48];
tdg q[48];
cx q[49], q[48];
t q[48];
cx q[49], q[4];
tdg q[4];
cx q[49], q[4];
t q[49];
t q[4];
h q[4];
cx q[67], q[16];
s q[30];
s q[25];
cx q[31], q[6];
h q[3];
cx q[3], q[44];
tdg q[44];
cx q[43], q[44];
t q[44];
cx q[3], q[44];
tdg q[44];
cx q[43], q[44];
t q[44];
cx q[43], q[3];
tdg q[3];
cx q[43], q[3];
t q[43];
t q[3];
h q[3];
h q[2];
s q[60];
h q[36];
cx q[36], q[43];
tdg q[43];
cx q[52], q[43];
t q[43];
cx q[36], q[43];
tdg q[43];
cx q[52], q[43];
t q[43];
cx q[52], q[36];
tdg q[36];
cx q[52], q[36];
t q[52];
t q[36];
h q[36];
h q[2];
s q[72];
h q[62];
h q[77];
cx q[77], q[10];
tdg q[10];
cx q[0], q[10];
t q[10];
cx q[77], q[10];
tdg q[10];
cx q[0], q[10];
t q[10];
cx q[0], q[77];
tdg q[77];
cx q[0], q[77];
t q[0];
t q[77];
h q[77];
h q[22];
cx q[22], q[32];
tdg q[32];
cx q[63], q[32];
t q[32];
cx q[22], q[32];
tdg q[32];
cx q[63], q[32];
t q[32];
cx q[63], q[22];
tdg q[22];
cx q[63], q[22];
t q[63];
t q[22];
h q[22];
t q[27];
s q[4];
s q[74];
s q[45];
cx q[21], q[75];
s q[76];
h q[39];
h q[64];
cx q[2], q[75];
s q[64];
t q[23];
t q[42];
h q[57];
cx q[57], q[26];
tdg q[26];
cx q[30], q[26];
t q[26];
cx q[57], q[26];
tdg q[26];
cx q[30], q[26];
t q[26];
cx q[30], q[57];
tdg q[57];
cx q[30], q[57];
t q[30];
t q[57];
h q[57];
h q[76];
cx q[76], q[30];
tdg q[30];
cx q[64], q[30];
t q[30];
cx q[76], q[30];
tdg q[30];
cx q[64], q[30];
t q[30];
cx q[64], q[76];
tdg q[76];
cx q[64], q[76];
t q[64];
t q[76];
h q[76];
h q[46];
cx q[46], q[45];
tdg q[45];
cx q[77], q[45];
t q[45];
cx q[46], q[45];
tdg q[45];
cx q[77], q[45];
t q[45];
cx q[77], q[46];
tdg q[46];
cx q[77], q[46];
t q[77];
t q[46];
h q[46];
cx q[61], q[63];
h q[28];
cx q[28], q[70];
tdg q[70];
cx q[5], q[70];
t q[70];
cx q[28], q[70];
tdg q[70];
cx q[5], q[70];
t q[70];
cx q[5], q[28];
tdg q[28];
cx q[5], q[28];
t q[5];
t q[28];
h q[28];
h q[38];
cx q[38], q[48];
tdg q[48];
cx q[40], q[48];
t q[48];
cx q[38], q[48];
tdg q[48];
cx q[40], q[48];
t q[48];
cx q[40], q[38];
tdg q[38];
cx q[40], q[38];
t q[40];
t q[38];
h q[38];
cx q[44], q[72];
h q[43];
cx q[18], q[72];
h q[46];
cx q[46], q[3];
tdg q[3];
cx q[59], q[3];
t q[3];
cx q[46], q[3];
tdg q[3];
cx q[59], q[3];
t q[3];
cx q[59], q[46];
tdg q[46];
cx q[59], q[46];
t q[59];
t q[46];
h q[46];
h q[1];
h q[78];
cx q[78], q[25];
tdg q[25];
cx q[23], q[25];
t q[25];
cx q[78], q[25];
tdg q[25];
cx q[23], q[25];
t q[25];
cx q[23], q[78];
tdg q[78];
cx q[23], q[78];
t q[23];
t q[78];
h q[78];
h q[49];
s q[8];
s q[23];
h q[10];
cx q[10], q[5];
tdg q[5];
cx q[38], q[5];
t q[5];
cx q[10], q[5];
tdg q[5];
cx q[38], q[5];
t q[5];
cx q[38], q[10];
tdg q[10];
cx q[38], q[10];
t q[38];
t q[10];
h q[10];
cx q[28], q[20];
h q[47];
cx q[56], q[61];
h q[59];
h q[4];
cx q[5], q[30];
h q[63];
cx q[63], q[20];
tdg q[20];
cx q[53], q[20];
t q[20];
cx q[63], q[20];
tdg q[20];
cx q[53], q[20];
t q[20];
cx q[53], q[63];
tdg q[63];
cx q[53], q[63];
t q[53];
t q[63];
h q[63];
t q[55];
cx q[4], q[27];
cx q[50], q[53];
h q[62];
cx q[62], q[31];
tdg q[31];
cx q[63], q[31];
t q[31];
cx q[62], q[31];
tdg q[31];
cx q[63], q[31];
t q[31];
cx q[63], q[62];
tdg q[62];
cx q[63], q[62];
t q[63];
t q[62];
h q[62];
s q[56];
s q[49];
h q[24];
h q[14];
cx q[14], q[37];
tdg q[37];
cx q[24], q[37];
t q[37];
cx q[14], q[37];
tdg q[37];
cx q[24], q[37];
t q[37];
cx q[24], q[14];
tdg q[14];
cx q[24], q[14];
t q[24];
t q[14];
h q[14];
h q[73];
s q[52];
h q[52];
cx q[52], q[44];
tdg q[44];
cx q[59], q[44];
t q[44];
cx q[52], q[44];
tdg q[44];
cx q[59], q[44];
t q[44];
cx q[59], q[52];
tdg q[52];
cx q[59], q[52];
t q[59];
t q[52];
h q[52];
t q[16];
t q[71];
h q[76];
h q[16];
cx q[29], q[24];
h q[23];
h q[15];
t q[0];
h q[60];
cx q[60], q[65];
tdg q[65];
cx q[56], q[65];
t q[65];
cx q[60], q[65];
tdg q[65];
cx q[56], q[65];
t q[65];
cx q[56], q[60];
tdg q[60];
cx q[56], q[60];
t q[56];
t q[60];
h q[60];
h q[49];
s q[19];
cx q[10], q[7];
s q[70];
t q[61];
s q[11];
h q[42];
t q[38];
cx q[20], q[47];
t q[45];
cx q[10], q[3];
t q[66];
h q[53];
cx q[53], q[22];
tdg q[22];
cx q[1], q[22];
t q[22];
cx q[53], q[22];
tdg q[22];
cx q[1], q[22];
t q[22];
cx q[1], q[53];
tdg q[53];
cx q[1], q[53];
t q[1];
t q[53];
h q[53];
s q[60];
cx q[46], q[44];
h q[47];
h q[73];
cx q[73], q[39];
tdg q[39];
cx q[27], q[39];
t q[39];
cx q[73], q[39];
tdg q[39];
cx q[27], q[39];
t q[39];
cx q[27], q[73];
tdg q[73];
cx q[27], q[73];
t q[27];
t q[73];
h q[73];
t q[0];
cx q[1], q[71];
h q[39];
h q[75];
h q[37];
s q[75];
h q[17];
h q[21];
t q[8];
s q[22];
cx q[31], q[22];
h q[25];
s q[70];
cx q[8], q[23];
s q[17];
cx q[59], q[26];
h q[76];
cx q[28], q[76];
t q[31];
t q[24];
h q[19];
t q[23];
h q[43];
cx q[43], q[68];
tdg q[68];
cx q[7], q[68];
t q[68];
cx q[43], q[68];
tdg q[68];
cx q[7], q[68];
t q[68];
cx q[7], q[43];
tdg q[43];
cx q[7], q[43];
t q[7];
t q[43];
h q[43];
h q[70];
s q[32];
s q[61];
t q[71];
t q[31];
h q[74];
cx q[74], q[54];
tdg q[54];
cx q[25], q[54];
t q[54];
cx q[74], q[54];
tdg q[54];
cx q[25], q[54];
t q[54];
cx q[25], q[74];
tdg q[74];
cx q[25], q[74];
t q[25];
t q[74];
h q[74];
h q[34];
t q[72];
h q[3];
cx q[3], q[75];
tdg q[75];
cx q[52], q[75];
t q[75];
cx q[3], q[75];
tdg q[75];
cx q[52], q[75];
t q[75];
cx q[52], q[3];
tdg q[3];
cx q[52], q[3];
t q[52];
t q[3];
h q[3];
t q[35];
cx q[73], q[48];
h q[77];
cx q[77], q[41];
tdg q[41];
cx q[49], q[41];
t q[41];
cx q[77], q[41];
tdg q[41];
cx q[49], q[41];
t q[41];
cx q[49], q[77];
tdg q[77];
cx q[49], q[77];
t q[49];
t q[77];
h q[77];
h q[11];
cx q[11], q[33];
tdg q[33];
cx q[12], q[33];
t q[33];
cx q[11], q[33];
tdg q[33];
cx q[12], q[33];
t q[33];
cx q[12], q[11];
tdg q[11];
cx q[12], q[11];
t q[12];
t q[11];
h q[11];
t q[64];
h q[50];
h q[68];
s q[39];
t q[23];
h q[42];
h q[24];
cx q[53], q[5];
h q[57];
h q[68];
h q[37];
cx q[37], q[34];
tdg q[34];
cx q[69], q[34];
t q[34];
cx q[37], q[34];
tdg q[34];
cx q[69], q[34];
t q[34];
cx q[69], q[37];
tdg q[37];
cx q[69], q[37];
t q[69];
t q[37];
h q[37];
cx q[41], q[42];
h q[75];
cx q[75], q[54];
tdg q[54];
cx q[79], q[54];
t q[54];
cx q[75], q[54];
tdg q[54];
cx q[79], q[54];
t q[54];
cx q[79], q[75];
tdg q[75];
cx q[79], q[75];
t q[79];
t q[75];
h q[75];
t q[56];
h q[62];
cx q[11], q[48];
s q[72];
t q[44];
