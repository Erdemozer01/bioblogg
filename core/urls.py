from django.contrib import admin
from django.urls import path, include
from django.views import generic
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from blog.views import ArchiveView

urlpatterns = [
                  path('admin/', admin.site.urls),
                  path('accounts/', include('django.contrib.auth.urls')),
                  path('laboratory/bioinformatic/', include('django_plotly_dash.urls')),
                  path('ckeditor/', include('ckeditor_uploader.urls')),
                  path('blog/', include("blog.urls", namespace="blog")),
                  path('accounts/', include("accounts.urls", namespace="accounts")),
                  path('dashboard/', include("dashboard.urls", namespace="dashboard")),
                  path('laboratory/bioinformatic/', include("bioinformatic.urls", namespace="bioinformatic")),
                  path('archive/<int:year>/<int:month>/', ArchiveView.as_view(), name="archive"),
                  path('', generic.TemplateView.as_view(template_name="cover.html"), name="anasayfa"),
              ] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

urlpatterns += staticfiles_urlpatterns()
