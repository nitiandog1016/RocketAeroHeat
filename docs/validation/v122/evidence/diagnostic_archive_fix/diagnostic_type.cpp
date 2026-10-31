#include "AbstractState.h"
#include "CoolProp.h"
#include "Exceptions.h"
#include <iostream>
#include <memory>
int main() {
  std::unique_ptr<CoolProp::AbstractState> s(CoolProp::AbstractState::factory("HEOS", "Air"));
  try { s->update(CoolProp::DmassP_INPUTS,5000.0/(287.0*4000.0),5000.0); }
  catch (const CoolProp::OutOfRangeError& e) { std::cout << "typed-OutOfRangeError " << e.what(); return 0; }
  catch (const std::exception& e) { std::cout << "other-exception " << e.what(); return 0; }
  std::cout << "unexpected-accept"; return 1;
}
