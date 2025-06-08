window.dash_clientside = Object.assign({}, window.dash_clientside, {
    molecule_viewer: {
        setupViewer: function (molData) {
            if (!molData) return;

            $('#viewer-container').empty();

            let $3dc = $('#viewer-container');
            let config = { backgroundColor: 'white' };
            let viewer = $3Dmol.createViewer($3dc, config);

            viewer.addModel(molData, "xyz");
            viewer.setStyle({ stick: { colorscheme: 'Jmol' } });
            viewer.zoomTo();
            viewer.render();
            viewer.zoom(1.2, 0);

            return '';
        }
    }
});
