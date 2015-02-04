import QtQuick 1.0

Rectangle {
	id: page
	width: 500; height: 200
	color: "lightgray"

	Text{
		id:helloText
		text: "Hello World!"
		y: 30
		anchors.horizontalCenter: page.horizontalCenter
		font.pointSize: 24; font.bold: true

		 MouseArea { id: mouseArea; anchors.fill: parent }

         states: State {
             name: "down"; when: mouseArea.pressed == true
             PropertyChanges { target: helloText; y: 160; rotation: 180; color: "red" }
         }

         transitions: Transition {
             from: ""; to: "down"; reversible: true
             ParallelAnimation {
                 NumberAnimation { properties: "y,rotation"; duration: 500; easing.type: Easing.InOutQuad }
                 ColorAnimation { duration: 500 }
             }
         }
	}

	Rectangle{
		id: cc
		width: 0
		height: 0
		x: 20
		y: 20
		color: "green"
		radius: width*0.5
		MouseArea {id: mm; anchors.fill: parent}

		states: State {
             name: "cc"; when: mouseArea.pressed == true
             PropertyChanges { target: cc; rotation: 180; x: 15; y: 15; height: 10; width: 10; color: "red" }
         }

         transitions: Transition {
             from: ""; to: "cc"; reversible: true
             ParallelAnimation {
                 NumberAnimation { properties: "x,y,rotation,height,width"; duration: 500; easing.type: Easing.InOutQuad }
                 ColorAnimation { duration: 500 }
             }
         }
	}

	Grid {
		id: colorPicker
		x: 4; anchors.bottom: page.bottom; anchors.bottomMargin: 4
		rows: 2; columns: 3; spacing: 3

		Cell { cellColor: "red"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "green"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "blue"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "yellow"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "steelblue"; onClicked: helloText.color = cellColor }
		Cell { cellColor: "black"; onClicked: helloText.color = cellColor }
	}
}