"""Page routes for the Enzyme Viewer web UI."""

from flask import render_template


def index():
    return render_template("index.html")


def favicon():
    return ("", 204)


def motif_library():
    return render_template("motif_library.html")


def test_nanozyme():
    return render_template("test_nanozyme.html")


def motif_view():
    return render_template("motif_view.html")


def nanozyme_design():
    return render_template("nanozyme_design.html")


def nanozyme_activity_validation():
    return render_template("activity_validation.html")


def register_page_routes(app) -> None:
    """Register UI page routes while preserving the original endpoint names."""
    app.add_url_rule("/", endpoint="index", view_func=index)
    app.add_url_rule("/favicon.ico", endpoint="favicon", view_func=favicon)
    app.add_url_rule("/motif_library", endpoint="motif_library", view_func=motif_library)
    app.add_url_rule("/test_nanozyme", endpoint="test_nanozyme", view_func=test_nanozyme)
    app.add_url_rule("/motif_view", endpoint="motif_view", view_func=motif_view)
    app.add_url_rule("/nanozyme_design", endpoint="nanozyme_design", view_func=nanozyme_design)
    app.add_url_rule(
        "/nanozyme_activity_validation",
        endpoint="nanozyme_activity_validation",
        view_func=nanozyme_activity_validation,
    )
