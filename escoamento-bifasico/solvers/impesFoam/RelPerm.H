// #include <cmath>

void printField(volScalarField field){
    forAll(field, i){
        Pout << field[i] << " ";
    }
    Pout << endl;
}

volScalarField correct_krw(const volScalarField Sw){
    // parameters 
    // double mu_w = 1.05e-3;
    double Swc = 0.1;
    double Sor = 0.15;
    double krwMax = 0.5;
    double a = 8;

    volScalarField Se = (Sw-Swc)/(1.0-Sw-Sor);

    // printField(Se);

    // forAll(Se,i){
    //     Se[i] = krw0 * std::pow(Se[i], nw);
    //     Pout << Se[i] << endl;
    // }

    return krwMax * Foam::pow(Se, a); // Sw;
}

volScalarField correct_kro(const volScalarField Sw){
    // parameters 
    // double mu_o = 1.05e-3;
    double Swc = 0.1;
    double Sor = 0.15;
    double kroMax = 1.0;
    double b = 3;

    volScalarField Se = (Sw-Swc)/(1.0-Swi-Sor);

    // printField(Se);
    
    // forAll(Se,i){
    //     Se[i] = krg0 * std::pow(1.0-Se[i], ng);
    // }

    return kroMax * Foam::pow(scalar(1.0)-Se, b); // 1.0-Sw;
}

volScalarField correct_fw(const volScalarField krw, const volScalarField kro){
    double mu_w = 1.0e-03;
    double mu_o = 20.0*mu_w;

    return (krw/mu_w)/(kro/mu_o + krw/mu_w);
    // return (1.0) / (1.0 + (mu_w*kro)/(mu_o*krw));
}

