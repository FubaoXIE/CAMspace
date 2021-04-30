#include "CAMGUIHelper.h"

#include "Entry.h"


void callback() { 

	ImGui::PushItemWidth(100);

	Basic_entry();
	my_entry();
	//make a difference 
	ImGui::SameLine();
	// ImGui::InputInt("source vertex", &iVertexSource);
	ImGui::PopItemWidth();
}

int main(int argc, char** argv) {

  // Options
  polyscope::options::programName = "CAMspace";
  polyscope::options::printPrefix = "[CAMspace] ";
  polyscope::options::autocenterStructures = false;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;

  // Initialize polyscope
  polyscope::init();

  // Add the callback
  polyscope::state::userCallback = callback;

  // Show the gui
  polyscope::show();

  return EXIT_SUCCESS;
}