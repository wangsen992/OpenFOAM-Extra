
#include "wrfFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFvPatchField<Type>(p, iF, fld)
{}


template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    fieldName_(dict.lookup<word>("field"))
{
    if (dict.found("value"))
    {
        Field<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        FatalIOErrorInFunction
        (
            dict
        )   << "Essential entry 'value' missing"
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const wrfFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper, false) // Don't map
{
    // Evaluate since value not mapped
    this->evaluate();
}


template<class Type>
Foam::wrfFixedValueFvPatchField<Type>::wrfFixedValueFvPatchField
(
    const wrfFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::wrfFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Temporary fix, doesn't support restarting simulations
    if(this->patch().boundaryMesh().mesh().time().timeIndex() > 0)
    {
      const Time& runTime = this->patch().boundaryMesh().mesh().time();
      WRF& wrf_ = runTime.template lookupObjectRef<WRF>("WRF");
    Info << "[fvPatchField] update coeffs for patch " << this->patch().name() << " for field " << fieldName_ << endl;
    const scalar t = this->db().time().timeOutputValue();
    typedef GeometricField<Type, fvPatchField, volMesh> psiType;
    psiType& psi
    (
      runTime.lookupObjectRef<psiType>(IOobject::groupName(fieldName_, "proj"))
    );
    // const pointField& Cf(this->patch().Cf());
    // Info << psi.boundaryField()[this->patch().index()].patchInternalField() << endl;

    this->operator==(psi.boundaryField()[this->patch().index()].patchInternalField());
    Info << average(*this) << endl;
    fixedValueFvPatchField<Type>::updateCoeffs();
    }
    else
    {
      return;
    }
    

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::wrfFixedValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    writeEntry(os, "field", fieldName_);
    writeEntry(os, "value", *this);
}


// ************************************************************************* //
